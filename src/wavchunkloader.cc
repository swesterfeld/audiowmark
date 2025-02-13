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

WavChunkLoader::WavChunkLoader (const std::string& filename, double chunk_size, double rate)
{
  m_filename = filename;
  m_chunk_size = chunk_size;
  m_rate = rate;
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

  m_wav_data = WavData ({}, m_in_stream->n_channels(), m_rate, m_in_stream->bit_depth());

  /* initialize resampler if input sample rate != watermark rate */
  m_input_rate = m_in_stream->sample_rate();
  if (m_input_rate != m_rate)
    m_resampler.reset (ResamplerImpl::create (m_in_stream->n_channels(), m_in_stream->sample_rate(), m_rate));

  /* samples2 buffer: a bit larger than the maximum length ClipDecoder uses (* speed_factor) */
  double speed_factor = 1.3;
  double block_seconds = (mark_sync_frame_count() + mark_data_frame_count()) * Params::frame_size / double (Params::mark_sample_rate);
  m_samples2_max_size = lrint (clip_decoder_max_blocks() * block_seconds * speed_factor * m_wav_data.sample_rate()) * m_wav_data.n_channels();
  m_samples2.reserve (m_samples2_max_size);

  /* samples1 buffer: chunk_size minutes + 5 minutes (to merge with samples2 buffer) */
  m_wav_data_max_size = m_samples2_max_size + lrint (m_chunk_size * 60 * m_wav_data.sample_rate()) * m_wav_data.n_channels();

  /* overlap size: BlockDecoder needs an overlap of 1 AB block (* speed factor) */
  m_n_overlap_samples = lrint (2 * block_seconds * speed_factor * m_wav_data.sample_rate()) * m_wav_data.n_channels();

  if (m_in_stream->n_frames() != AudioInputStream::N_FRAMES_UNKNOWN)
    {
      size_t n_reserve_frames = m_in_stream->n_frames() * double (m_rate) / m_input_rate;

      /* since we possibly resample the input signal, we expect a slight difference between
       * predicted and actual space requirements, so we reserve a bit more than predicted
       */
      n_reserve_frames *= 1.001;
      n_reserve_frames += 100;

      m_wav_data.mutable_samples().reserve (std::min (m_wav_data_max_size, n_reserve_frames * m_wav_data.n_channels()));
    }
  else
    {
      // unknown size: only reserve 5 minutes (not chunk_size) in order to minimize memory usage for short files
      m_wav_data.mutable_samples().reserve (m_samples2_max_size);
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

  if (m_wav_data.n_values()) // avoid division by zero for empty wav_data
    m_time_offset += m_wav_data.n_frames() / double (m_wav_data.sample_rate());

  vector<float>& ref_samples1 = m_wav_data.mutable_samples();

  if (!ref_samples1.empty()) /* second block or later */
    {
      /* overlap samples1 with last block */
      assert (ref_samples1.size() >= m_n_overlap_samples);
      ref_samples1.erase (ref_samples1.begin(), ref_samples1.end() - m_n_overlap_samples);

      m_time_offset -= (m_n_overlap_samples / m_wav_data.n_channels()) / double (m_wav_data.sample_rate());

      /* append samples2 to samples1 */
      ref_samples1.insert (ref_samples1.end(), m_samples2.begin(), m_samples2.end());
      m_samples2.clear();
    }

  bool eof = false;
  Error err = refill (ref_samples1, m_wav_data_max_size - m_samples2_max_size, m_wav_data_max_size, &eof);
  if (err)
    {
      m_state = State::ERROR;
      return err;
    }

  if (!eof)
    {
      err = refill (m_samples2, m_samples2_max_size, m_samples2_max_size, &eof);
      if (err)
        {
          m_state = State::ERROR;
          return err;
        }
    }

  if (eof)
    {
      /* for long files:
       *  - ensure that last chunk is large enough -> avoid ClipDecoder
      */
      ref_samples1.insert (ref_samples1.end(), m_samples2.begin(), m_samples2.end());
      m_samples2.clear();

      if (ref_samples1.size())
        m_state = State::LAST_CHUNK;
      else
        m_state = State::DONE;
    }

  if (Params::test_truncate)
    {
      const size_t want_n_samples = m_wav_data.sample_rate() * m_wav_data.n_channels() * Params::test_truncate;
      if (want_n_samples > m_wav_data_max_size - m_samples2_max_size)
        return Error ("test truncate must be less than chunk size");

      if (want_n_samples < ref_samples1.size())
        ref_samples1.resize (want_n_samples);

      if (ref_samples1.size())
        m_state = State::LAST_CHUNK;
      else
        m_state = State::DONE;
    }

  printf ("chunk size: %f minutes, cap %f minutes and %f minutes\n", ref_samples1.size() / m_in_stream->n_channels() / (60.0 * m_rate),
                                                                     ref_samples1.capacity() / m_in_stream->n_channels() / (60.0 * m_rate),
                                                                     m_samples2.capacity() / m_in_stream->n_channels() / (60.0 * m_rate));
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
WavChunkLoader::refill (std::vector<float>& samples, size_t values, size_t max_size, bool *eof)
{
  *eof = false;

  constexpr size_t block_size = 4096;

  vector<float> buffer;
  while (samples.size() < values)
    {
      if (m_resampler)
        {
          if (m_resampler->can_read_frames() < block_size && !m_resampler_in_eof)
            {
              Error err = m_in_stream->read_frames (buffer, block_size * double (m_input_rate) / m_rate);
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

          buffer = m_resampler->read_frames (std::min<size_t> (m_resampler->can_read_frames(), (values - samples.size()) / m_wav_data.n_channels()));
        }
      else
        {
          Error err = m_in_stream->read_frames (buffer, std::min<size_t> (block_size, (values - samples.size()) / m_wav_data.n_channels()));
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
