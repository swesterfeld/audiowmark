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

#ifndef AUDIOWMARK_WAV_PIPE_INPUT_STREAM_HH
#define AUDIOWMARK_WAV_PIPE_INPUT_STREAM_HH

#include <string>
#include <memory>

#include <sndfile.h>

#include "audiostream.hh"
#include "rawinputstream.hh"

class WavPipeInputStream : public AudioInputStream
{
  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  State       m_state = State::NEW;
  RawFormat   m_format;
  FILE       *m_input_file = nullptr;
  bool        m_close_file = false;

  std::vector<unsigned char>    m_input_bytes;
  std::unique_ptr<RawConverter> m_raw_converter;

  Error read_error (const std::string& message);

public:
  ~WavPipeInputStream();

  Error   open (const std::string& filename);
  Error   read_frames (std::vector<float>& samples, size_t count) override;
  void    close();

  int     bit_depth() const override;
  int     sample_rate() const override;
  size_t  n_frames() const override;
  int     n_channels() const override;
};

#endif /* AUDIOWMARK_WAV_PIPE_INPUT_STREAM_HH */
