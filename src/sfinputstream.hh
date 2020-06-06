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

#ifndef AUDIOWMARK_SF_INPUT_STREAM_HH
#define AUDIOWMARK_SF_INPUT_STREAM_HH

#include <string>

#include <sndfile.h>

#include "audiostream.hh"

class SFInputStream : public AudioInputStream
{
  SNDFILE    *m_sndfile = nullptr;
  int         m_n_channels = 0;
  int         m_n_values = 0;
  int         m_bit_depth = 0;
  int         m_sample_rate = 0;

  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  State       m_state = State::NEW;

public:
  ~SFInputStream();

  Error               open (const std::string& filename);
  Error               read_frames (std::vector<float>& samples, size_t count) override;
  void                close();

  int
  n_channels() const override
  {
    return m_n_channels;
  }
  int sample_rate() const override;
  int bit_depth() const override;
  size_t
  n_values() const
  {
    return m_n_values;
  }
  size_t
  n_frames() const override
  {
    return m_n_values / m_n_channels;
  }
};

#endif /* AUDIOWMARK_SF_INPUT_STREAM_HH */
