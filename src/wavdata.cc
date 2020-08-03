/*
 * Copyright (C) 2018-2020 Stefan Westerfeld
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

#include "wavdata.hh"
#include "utils.hh"
#include "audiostream.hh"
#include "sfinputstream.hh"
#include "sfoutputstream.hh"
#include "mp3inputstream.hh"

#include <memory>
#include <math.h>

using std::string;
using std::vector;

WavData::WavData()
{
}

WavData::WavData (const vector<float>& samples, int n_channels, int sample_rate, int bit_depth)
{
  m_samples     = samples;
  m_n_channels  = n_channels;
  m_sample_rate = sample_rate;
  m_bit_depth   = bit_depth;
}

Error
WavData::load (const string& filename)
{
  Error err;

  std::unique_ptr<AudioInputStream> in_stream = AudioInputStream::create (filename, err);
  if (err)
    return err;

  return load (in_stream.get());
}

Error
WavData::load (AudioInputStream *in_stream)
{
  m_samples.clear(); // get rid of old contents

  vector<float> m_buffer;
  while (true)
    {
      Error err = in_stream->read_frames (m_buffer, 1024);
      if (err)
        return err;

      if (!m_buffer.size())
        {
          /* reached eof */
          break;
        }
      m_samples.insert (m_samples.end(), m_buffer.begin(), m_buffer.end());
    }
  m_sample_rate = in_stream->sample_rate();
  m_n_channels  = in_stream->n_channels();
  m_bit_depth   = in_stream->bit_depth();

  return Error::Code::NONE;
}

Error
WavData::save (const string& filename) const
{
  std::unique_ptr<AudioOutputStream> out_stream;
  Error err;

  out_stream = AudioOutputStream::create (filename, m_n_channels, m_sample_rate, m_bit_depth, m_samples.size() / m_n_channels, err);
  if (err)
    return err;

  err = out_stream->write_frames (m_samples);
  if (err)
    return err;

  err = out_stream->close();
  return err;
}

int
WavData::sample_rate() const
{
  return m_sample_rate;
}

int
WavData::bit_depth() const
{
  return m_bit_depth;
}

void
WavData::set_samples (const vector<float>& samples)
{
  m_samples = samples;
}
