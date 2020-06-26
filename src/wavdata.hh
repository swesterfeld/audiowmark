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

#ifndef AUDIOWMARK_WAV_DATA_HH
#define AUDIOWMARK_WAV_DATA_HH

#include <string>
#include <vector>

#include "utils.hh"
#include "audiostream.hh"

class WavData
{
  std::vector<float> m_samples;
  int                m_sample_rate = 0;
  int                m_n_channels  = 0;
  int                m_bit_depth   = 0;

public:
  WavData();
  WavData (const std::vector<float>& samples, int n_channels, int sample_rate, int bit_depth);

  Error load (AudioInputStream *in_stream);
  Error load (const std::string& filename);
  Error save (const std::string& filename) const;

  int                         sample_rate() const;
  int                         bit_depth() const;

  int
  n_channels() const
  {
    return m_n_channels;
  }
  size_t
  n_values() const
  {
    return m_samples.size();
  }
  size_t
  n_frames() const
  {
    return m_samples.size() / m_n_channels;
  }
  const std::vector<float>&
  samples() const
  {
    return m_samples;
  }

  void set_samples (const std::vector<float>& samples);
};

#endif /* AUDIOWMARK_WAV_DATA_HH */
