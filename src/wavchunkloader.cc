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
  if (!m_in_stream)
    {
      Error err;

      m_in_stream = AudioInputStream::create (m_filename, err);
      if (err)
        return err;
    }

  if (m_wav_data.n_values()) // avoid division by zero for empty wav_data
    m_time_offset += m_wav_data.n_frames() / double (m_wav_data.sample_rate());

  size_t overlap = m_wav_data.sample_rate() * m_wav_data.n_channels() * 60 * 5;
  printf ("overlap=%zd\n", overlap);
  if (m_samples1.size() > overlap)
    {
      m_time_offset -= overlap / m_wav_data.sample_rate() / m_wav_data.n_channels();
      m_samples1.erase (m_samples1.begin(), m_samples1.end() - overlap);
    }
  else
    {
      m_samples1.clear();
    }
  m_samples1.insert (m_samples1.end(), m_samples2.begin(), m_samples2.end());
  m_samples2.clear();

  bool eof = !refill (m_samples1, m_chunk_size);
  if (!eof)
    eof = !refill (m_samples2, 5);

  if (eof)
    {
      /* for long files:
       *  - ensure that last chunk is large enough -> avoid ClipDecoder
      */
      m_samples1.insert (m_samples1.end(), m_samples2.begin(), m_samples2.end());
      m_samples2.clear();
    }

  m_wav_data = WavData (m_samples1, m_in_stream->n_channels(), m_in_stream->sample_rate(), m_in_stream->bit_depth());
  printf ("chunk size: %f minutes\n", m_samples1.size() / m_in_stream->n_channels() / (60.0 * m_in_stream->sample_rate()));

  return Error::Code::NONE;
}

bool
WavChunkLoader::refill (std::vector<float>& samples, double time)
{
  auto to_minutes = [&] (size_t n_samples) {
    return (n_samples / m_in_stream->n_channels()) / (60.0 * m_in_stream->sample_rate());
  };

  vector<float> m_buffer;
  while (to_minutes (samples.size()) < time)
    {
      Error err = m_in_stream->read_frames (m_buffer, 1024);
      if (err)
        return err;

      if (!m_buffer.size())
        {
          /* reached eof */
          return false;
        }
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
