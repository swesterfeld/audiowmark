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

#include "rawinputstream.hh"
#include "rawconverter.hh"

#include <assert.h>
#include <string.h>
#include <errno.h>

using std::string;
using std::vector;

RawFormat::RawFormat()
{
}

RawFormat::RawFormat (int n_channels, int sample_rate, int bit_depth) :
  m_n_channels (n_channels),
  m_sample_rate (sample_rate),
  m_bit_depth (bit_depth)
{
}

void
RawFormat::set_channels (int channels)
{
  m_n_channels = channels;
}

void
RawFormat::set_sample_rate (int rate)
{
  m_sample_rate = rate;
}

void
RawFormat::set_bit_depth (int bits)
{
  m_bit_depth = bits;
}

void
RawFormat::set_endian (Endian endian)
{
  m_endian = endian;
}

void
RawFormat::set_encoding (Encoding encoding)
{
  m_encoding = encoding;
}

RawInputStream::~RawInputStream()
{
  close();
}

Error
RawInputStream::open (const string& filename, const RawFormat& format)
{
  assert (m_state == State::NEW);

  if (!format.n_channels())
    return Error ("RawInputStream: input format: missing number of channels");
  if (!format.bit_depth())
    return Error ("RawInputStream: input format: missing bit depth");
  if (!format.sample_rate())
    return Error ("RawInputStream: input format: missing sample rate");

  Error err = Error::Code::NONE;
  m_raw_converter.reset (RawConverter::create (format, err));
  if (err)
    return err;

  if (filename == "-")
    {
      m_input_file = stdin;
      m_close_file = false;
    }
  else
    {
      m_input_file = fopen (filename.c_str(), "r");
      if (!m_input_file)
        return Error (strerror (errno));

      m_close_file = true;
    }

  m_format = format;
  m_state  = State::OPEN;
  return Error::Code::NONE;
}

int
RawInputStream::sample_rate() const
{
  return m_format.sample_rate();
}

int
RawInputStream::bit_depth() const
{
  return m_format.bit_depth();
}

size_t
RawInputStream::n_frames() const
{
  return N_FRAMES_UNKNOWN;
}

int
RawInputStream::n_channels() const
{
  return m_format.n_channels();
}

Error
RawInputStream::read_frames (vector<float>& samples, size_t count)
{
  assert (m_state == State::OPEN);

  const int n_channels   = m_format.n_channels();
  const int sample_width = m_format.bit_depth() / 8;

  vector<unsigned char> input_bytes (count * n_channels * sample_width);
  size_t r_count = fread (input_bytes.data(), n_channels * sample_width, count, m_input_file);
  if (ferror (m_input_file))
    return Error ("error reading sample data");

  input_bytes.resize (r_count * n_channels * sample_width);

  m_raw_converter->from_raw (input_bytes, samples);

  return Error::Code::NONE;
}

void
RawInputStream::close()
{
  if (m_state == State::OPEN)
    {
      if (m_close_file && m_input_file)
        {
          fclose (m_input_file);
          m_input_file = nullptr;
          m_close_file = false;
        }

      m_state = State::CLOSED;
    }
}
