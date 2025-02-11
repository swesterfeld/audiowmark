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

#include <math.h>
#include <assert.h>

using std::vector;

WavChunkLoader::WavChunkLoader (const std::string& filename, double chunk_size, double rate)
{
  m_filename = filename;
  m_chunk_size = chunk_size;
  m_rate = rate;
  m_time_offset = 0;
}

Error
WavChunkLoader::load_next_chunk()
{
  vector<float>& ref_samples1 = m_wav_data.mutable_samples();

  if (!m_in_stream)
    {
      Error err;

      m_in_stream = AudioInputStream::create (m_filename, err);
      if (err)
        return err;


      m_wav_data = WavData ({}, m_in_stream->n_channels(), m_in_stream->sample_rate(), m_in_stream->bit_depth());

      /* samples2 buffer: 5 minutes */
      m_samples2_max_size = lrint (5 * 60 * m_wav_data.sample_rate()) * m_wav_data.n_channels();
      m_samples2.reserve (m_samples2_max_size);

      /* samples1 buffer: chunk_size minutes + 5 minutes (to merge with samples2 buffer) */
      m_wav_data_max_size = m_samples2_max_size + lrint (m_chunk_size * 60 * m_wav_data.sample_rate()) * m_wav_data.n_channels();

      // only reserve 5 minutes (not chunk_size) in order to minimize memory usage for short files
      ref_samples1.reserve (m_samples2_max_size);
    }

  if (m_wav_data.n_values()) // avoid division by zero for empty wav_data
    m_time_offset += m_wav_data.n_frames() / double (m_wav_data.sample_rate());

  size_t overlap = m_wav_data.sample_rate() * m_wav_data.n_channels() * 60 * 5;
  printf ("overlap=%zd\n", overlap);
  if (ref_samples1.size() > overlap)
    {
      m_time_offset -= overlap / m_wav_data.sample_rate() / m_wav_data.n_channels();
      ref_samples1.erase (ref_samples1.begin(), ref_samples1.end() - overlap);
    }
  else
    {
      ref_samples1.clear();
    }
  ref_samples1.insert (ref_samples1.end(), m_samples2.begin(), m_samples2.end());
  m_samples2.clear();

  bool eof = !refill (ref_samples1, m_wav_data_max_size - m_samples2_max_size, m_wav_data_max_size);
  if (!eof)
    eof = !refill (m_samples2, m_samples2_max_size, m_samples2_max_size);

  if (eof)
    {
      /* for long files:
       *  - ensure that last chunk is large enough -> avoid ClipDecoder
      */
      ref_samples1.insert (ref_samples1.end(), m_samples2.begin(), m_samples2.end());
      m_samples2.clear();
    }

  printf ("chunk size: %f minutes, cap %f minutes and %f minutes\n", ref_samples1.size() / m_in_stream->n_channels() / (60.0 * m_in_stream->sample_rate()),
                                                                     ref_samples1.capacity() / m_in_stream->n_channels() / (60.0 * m_in_stream->sample_rate()),
                                                                     m_samples2.capacity() / m_in_stream->n_channels() / (60.0 * m_in_stream->sample_rate()));
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

bool
WavChunkLoader::refill (std::vector<float>& samples, size_t values, size_t max_size)
{
  vector<float> m_buffer;
  while (samples.size() < values)
    {
      Error err = m_in_stream->read_frames (m_buffer, std::min<size_t> (1024, (values - samples.size()) / m_wav_data.n_channels()));
      if (err)
        return err;

      if (!m_buffer.size())
        {
          /* reached eof */
          return false;
        }

      update_capacity (samples, samples.size() + m_buffer.size(), max_size);
      samples.insert (samples.end(), m_buffer.begin(), m_buffer.end());
    }
  return true;
}

bool
WavChunkLoader::done()
{
  return m_wav_data.n_frames() == 0;
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
