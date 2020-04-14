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

#ifndef AUDIOWMARK_RAW_INPUT_STREAM_HH
#define AUDIOWMARK_RAW_INPUT_STREAM_HH

#include <string>
#include <memory>

#include <sndfile.h>

#include "audiostream.hh"

class RawFormat
{
public:
  enum Endian {
    LITTLE,
    BIG
  };
  enum Encoding {
    SIGNED,
    UNSIGNED
  };
private:
  int       m_n_channels  = 2;
  int       m_sample_rate = 0;
  int       m_bit_depth   = 16;
  Endian    m_endian      = LITTLE;
  Encoding  m_encoding    = SIGNED;
public:
  RawFormat();
  RawFormat (int n_channels, int sample_rate, int bit_depth);

  int n_channels() const { return m_n_channels; }
  int sample_rate() const { return m_sample_rate; }
  int bit_depth() const { return m_bit_depth; }
  Endian endian() const { return m_endian; }
  Encoding encoding() const { return m_encoding; }

  void set_channels (int channels);
  void set_sample_rate (int rate);
  void set_bit_depth (int bits);
  void set_endian (Endian endian);
  void set_encoding (Encoding encoding);
};

class RawConverter;

class RawInputStream : public AudioInputStream
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

  std::unique_ptr<RawConverter> m_raw_converter;

public:
  ~RawInputStream();

  Error   open (const std::string& filename, const RawFormat& format);
  Error   read_frames (std::vector<float>& samples, size_t count) override;
  void    close();

  int     bit_depth() const override;
  int     sample_rate() const override;
  size_t  n_frames() const override;
  int     n_channels() const override;
};

#endif /* AUDIOWMARK_RAW_INPUT_STREAM_HH */
