/*
 * Copyright (C) 2025 Stefan Westerfeld
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "wavchunkloader.hh"
#include "wmcommon.hh"

#include <math.h>
#include <assert.h>

using std::vector;

/* WavChunkLoader reads the input file, resamples it to the watermark sample
 * rate and splits it into overlapping chunks. This works without knowing the
 * length of the input file before eof is reached.
 *
 * Example:
 *          <-------- m_wav_data --------->
 *          <-overlap->
 * Chunk 1: [ A A A A A A A A A A A A A A ]
 *
 * Chunk 2: [ A A A A A B B B B B B B B B ]
 *
 * Chunk 3: [ B B B B B C C ]
 *
 * Chunk 1: Fill m_wav_data by reading from the input stream (A).
 *
 * Chunk 2: Keep some overlap samples from the end of m_wav_data (A) and
 * fill the remaining space by reading from the input stream (B).
 *
 * Chunk 3: Keep some overlap samples (B) and try to fill the block from the
 * input stream (C). Since EOF was reached, we're done at this point.
 *
 * Overlap should be larger than one AB block to get all BlockDecoder results.
 */
WavChunkLoader::WavChunkLoader (const std::string& filename) :
  m_filename (filename)
{
}

Error
WavChunkLoader::open()
{
  assert (m_state == State::NEW);

  Error err;
  m_in_stream = AudioInputStream::create (m_filename, err);
  if (err)
    {
      m_state = State::ERROR;
      return err;
    }
  m_state = State::OPEN;

  m_wav_data = WavData ({}, m_in_stream->n_channels(), Params::mark_sample_rate, m_in_stream->bit_depth());

  /* initialize resampler if input sample rate != watermark rate */
  if (m_in_stream->sample_rate() != m_wav_data.sample_rate())
    m_resampler.reset (ResamplerImpl::create (m_in_stream->n_channels(), m_in_stream->sample_rate(), m_wav_data.sample_rate()));

  /* maximum length of the m_wav_data samples (chunk size) */
  m_wav_data_max_size = lrint (Params::get_chunk_size * 60 * m_wav_data.sample_rate()) * m_wav_data.n_channels();

  /* overlap size:
   *  - should be large enough for BlockDecoder overlap (1 AB block == 2 blocks)
   *  - take speed factor into account for speed detection
   */
  const int overlap_blocks = 2;
  const double speed_factor = 1.3;
  const double block_seconds = (mark_sync_frame_count() + mark_data_frame_count()) * Params::frame_size / double (Params::mark_sample_rate);
  m_n_overlap_samples = lrint (overlap_blocks * block_seconds * speed_factor * m_wav_data.sample_rate()) * m_wav_data.n_channels();

  if (m_in_stream->n_frames() != AudioInputStream::N_FRAMES_UNKNOWN)
    {
      size_t n_reserve_frames = m_in_stream->n_frames() * double (m_wav_data.sample_rate()) / m_in_stream->sample_rate();

      /* since we possibly resample the input signal, we expect a slight difference between
       * predicted and actual space requirements, so we reserve a bit more than predicted
       */
      n_reserve_frames *= 1.001;
      n_reserve_frames += 100;

      m_wav_data.mutable_samples().reserve (std::min (m_wav_data_max_size, n_reserve_frames * m_wav_data.n_channels()));
    }
  return Error::Code::NONE;
}

Error
WavChunkLoader::load_next_chunk()
{
  assert (m_state != State::ERROR);

  if (m_state == State::LAST_CHUNK)
    {
      m_state = State::DONE;
      return Error::Code::NONE;
    }

  if (m_state == State::NEW)
    {
      Error err = open();
      assert (m_state != State::NEW); // m_state should be State::OPEN or State::ERROR now

      if (err)
        return err;
    }

  vector<float>& ref_samples = m_wav_data.mutable_samples();
  if (!ref_samples.empty()) /* second block or later */
    {
      /* overlap samples with last block */
      assert (ref_samples.size() >= m_n_overlap_samples);

      m_time_offset += ((ref_samples.size() - m_n_overlap_samples) / m_wav_data.n_channels()) / double (m_wav_data.sample_rate());
      ref_samples.erase (ref_samples.begin(), ref_samples.end() - m_n_overlap_samples);
    }

  bool eof = false;
  Error err = refill (ref_samples, m_wav_data_max_size, &eof);
  if (err)
    {
      m_state = State::ERROR;
      return err;
    }

  if (eof)
    {
      if (ref_samples.size())
        m_state = State::LAST_CHUNK;
      else
        m_state = State::DONE;
    }

  if (Params::test_truncate)
    {
      const size_t want_n_samples = m_wav_data.sample_rate() * m_wav_data.n_channels() * Params::test_truncate;
      if (want_n_samples > m_wav_data_max_size)
        return Error ("test truncate must be less than chunk size");

      if (want_n_samples < ref_samples.size())
        ref_samples.resize (want_n_samples);

      if (ref_samples.size())
        m_state = State::LAST_CHUNK;
      else
        m_state = State::DONE;
    }

  return Error::Code::NONE;
}

void
WavChunkLoader::update_capacity (vector<float>& samples, size_t need_space, size_t max_size)
{
  assert (need_space <= max_size);

  if (samples.capacity() < need_space)
    {
      /* We double the capacity of the std::vector at each step. However, if
       * the new capacity exceeds 40% of the total maximum size, we use the
       * maximum total size instead. This approach reduces peak memory usage
       * compared to always doubling.
       */
      size_t cap = 8192;
      while (cap < need_space)
        cap *= 2;

      double new_percent_full = cap * 100. / max_size;
      if (new_percent_full > 40)
        cap = max_size;

      samples.reserve (cap);
    }

  assert (samples.capacity() >= need_space);
}

Error
WavChunkLoader::refill (std::vector<float>& samples, size_t max_size, bool *eof)
{
  *eof = false;

  constexpr size_t block_size = 4096;

  vector<float> buffer;
  while (samples.size() < max_size)
    {
      if (m_resampler)
        {
          if (m_resampler->can_read_frames() < block_size && !m_resampler_in_eof)
            {
              Error err = m_in_stream->read_frames (buffer, block_size * double (m_in_stream->sample_rate()) / m_wav_data.sample_rate());
              if (err)
                return err;

              m_resampler->write_frames (buffer);
              if (!buffer.size())
                {
                  /* input file reached eof */
                  m_resampler->write_trailing_frames();
                  m_resampler_in_eof = true;
                }
            }

          buffer = m_resampler->read_frames (std::min<size_t> (m_resampler->can_read_frames(), (max_size - samples.size()) / m_wav_data.n_channels()));
        }
      else
        {
          Error err = m_in_stream->read_frames (buffer, std::min<size_t> (block_size, (max_size - samples.size()) / m_wav_data.n_channels()));
          if (err)
            return err;
        }

      if (!buffer.size())
        {
          /* reached eof */
          *eof = true;
          return Error::Code::NONE;
        }

      update_capacity (samples, samples.size() + buffer.size(), max_size);
      samples.insert (samples.end(), buffer.begin(), buffer.end());
      m_n_total_samples += buffer.size();
    }
  return Error::Code::NONE;
}

bool
WavChunkLoader::done()
{
  return m_state == State::DONE;
}

double
WavChunkLoader::length()
{
  assert (m_state == State::DONE);

  return m_n_total_samples / double (m_wav_data.sample_rate() * m_wav_data.n_channels());
}

const WavData&
WavChunkLoader::wav_data()
{
  return m_wav_data;
}

double
WavChunkLoader::time_offset()
{
  return m_time_offset;
}
