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

#ifndef AUDIOWMARK_STDOUT_WAV_STREAM_HH
#define AUDIOWMARK_STDOUT_WAV_STREAM_HH

#include "audiostream.hh"
#include "rawconverter.hh"

#include <string>

class StdoutWavOutputStream : public AudioOutputStream
{
  int         m_bit_depth = 0;
  int         m_sample_rate = 0;
  int         m_n_channels = 0;
  size_t      m_close_padding = 0;

  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  State       m_state = State::NEW;

  std::unique_ptr<RawConverter> m_raw_converter;

public:
  ~StdoutWavOutputStream();

  Error open (int n_channels, int sample_rate, int bit_depth, size_t n_frames);
  Error write_frames (const std::vector<float>& frames) override;
  Error close() override;
  int  sample_rate() const override;
  int  bit_depth() const override;
  int  n_channels() const override;
};

#endif
