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

#ifndef AUDIOWMARK_SF_OUTPUT_STREAM_HH
#define AUDIOWMARK_SF_OUTPUT_STREAM_HH

#include <string>

#include <sndfile.h>

#include "audiostream.hh"
#include "sfinputstream.hh"

class SFOutputStream : public AudioOutputStream
{
public:
  enum class OutFormat { WAV, RF64, FLAC };

private:
  SFVirtualData m_virtual_data;

  SNDFILE    *m_sndfile = nullptr;
  int         m_bit_depth = 0;
  int         m_sample_rate = 0;
  int         m_n_channels = 0;

  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  State       m_state = State::NEW;

  Error open (std::function<SNDFILE* (SF_INFO *)> open_func, int n_channels, int sample_rate, int bit_depth, OutFormat out_format);
public:
  ~SFOutputStream();

  Error  open (const std::string& filename, int n_channels, int sample_rate, int bit_depth, OutFormat out_format = OutFormat::WAV);
  Error  open (std::vector<unsigned char> *data, int n_channels, int sample_rate, int bit_depth, OutFormat out_format = OutFormat::WAV);
  Error  write_frames (const std::vector<float>& frames) override;
  Error  close() override;
  int    bit_depth() const override;
  int    sample_rate() const override;
  int    n_channels() const override;
};

#endif /* AUDIOWMARK_SF_OUTPUT_STREAM_HH */
