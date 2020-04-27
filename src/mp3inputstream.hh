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

#ifndef AUDIOWMARK_MP3_INPUT_STREAM_HH
#define AUDIOWMARK_MP3_INPUT_STREAM_HH

#include <string>
#include <mpg123.h>

#include "audiostream.hh"

class MP3InputStream : public AudioInputStream
{
  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  int         m_n_values = 0;
  int         m_n_channels = 0;
  int         m_sample_rate = 0;
  size_t      m_frames_left = 0;
  bool        m_need_close = false;
  bool        m_eof = false;
  State       m_state = State::NEW;

  mpg123_handle     *m_handle = nullptr;
  std::vector<float> m_read_buffer;
public:
  ~MP3InputStream();

  Error   open (const std::string& filename);
  Error   read_frames (std::vector<float>& samples, size_t count) override;
  void    close();

  int     bit_depth() const override;
  int     sample_rate() const override;
  int     n_channels()  const override;
  size_t  n_frames() const override;

  static bool detect (const std::string& filename);
};

#endif /* AUDIOWMARK_MP3_INPUT_STREAM_HH */
