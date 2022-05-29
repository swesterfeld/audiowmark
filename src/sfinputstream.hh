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
#include <functional>

#include <sndfile.h>

#include "audiostream.hh"

/* to support virtual io read/write from/to memory */
struct SFVirtualData
{
  SFVirtualData();

  std::vector<unsigned char> *mem    = nullptr;
  sf_count_t                  offset = 0;
  SF_VIRTUAL_IO               io;
};

class SFInputStream : public AudioInputStream
{
private:
  SFVirtualData m_virtual_data;

  SNDFILE    *m_sndfile = nullptr;
  int         m_n_channels = 0;
  size_t      m_n_frames = 0;
  int         m_bit_depth = 0;
  int         m_sample_rate = 0;
  bool        m_read_float_data = false;
  bool        m_is_stdin = false;

  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  State       m_state = State::NEW;

  Error open (std::function<SNDFILE* (SF_INFO *)> open_func);
public:
  ~SFInputStream();

  Error               open (const std::string& filename);
  Error               open (const std::vector<unsigned char> *data);
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
  n_frames() const override
  {
    return m_n_frames;
  }
};

#endif /* AUDIOWMARK_SF_INPUT_STREAM_HH */
