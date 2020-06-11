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

#ifndef AUDIOWMARK_RAW_OUTPUT_STREAM_HH
#define AUDIOWMARK_RAW_OUTPUT_STREAM_HH

#include "rawinputstream.hh"
#include "rawconverter.hh"

#include <memory>

class RawOutputStream : public AudioOutputStream
{
  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  State       m_state = State::NEW;
  RawFormat   m_format;
  FILE       *m_output_file = nullptr;
  bool        m_close_file = false;

  std::unique_ptr<RawConverter> m_raw_converter;
public:
  ~RawOutputStream();

  int   bit_depth() const override;
  int   sample_rate() const override;
  int   n_channels()  const override;

  Error open (const std::string& filename, const RawFormat& format);
  Error write_frames (const std::vector<float>& frames) override;
  Error close() override;
};

#endif /* AUDIOWMARK_RAW_OUTPUT_STREAM_HH */
