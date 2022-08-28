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

#include "rawoutputstream.hh"

#include <assert.h>
#include <string.h>
#include <errno.h>

using std::string;
using std::vector;

RawOutputStream::~RawOutputStream()
{
  close();
}

Error
RawOutputStream::open (const string& filename, const RawFormat& format)
{
  assert (m_state == State::NEW);

  if (!format.n_channels())
    return Error ("RawOutputStream: output format: missing number of channels");
  if (!format.bit_depth())
    return Error ("RawOutputStream: output format: missing bit depth");
  if (!format.sample_rate())
    return Error ("RawOutputStream: output format: missing sample rate");

  Error err = Error::Code::NONE;
  m_raw_converter.reset (RawConverter::create (format, err));
  if (err)
    return err;

  if (filename == "-")
    {
      m_output_file = stdout;
      m_close_file = false;
    }
  else
    {
      m_output_file = fopen (filename.c_str(), "w");
      if (!m_output_file)
        return Error (strerror (errno));

      m_close_file = true;
    }

  m_format = format;
  m_state  = State::OPEN;
  return Error::Code::NONE;
}

int
RawOutputStream::sample_rate() const
{
  return m_format.sample_rate();
}

int
RawOutputStream::bit_depth() const
{
  return m_format.bit_depth();
}

int
RawOutputStream::n_channels() const
{
  return m_format.n_channels();
}

Error
RawOutputStream::write_frames (const vector<float>& samples)
{
  assert (m_state == State::OPEN);

  if (samples.empty())
    return Error::Code::NONE;

  vector<unsigned char> bytes;
  m_raw_converter->to_raw (samples, bytes);

  fwrite (&bytes[0], 1, bytes.size(), m_output_file);
  if (ferror (m_output_file))
    return Error ("write sample data failed");

  return Error::Code::NONE;
}

Error
RawOutputStream::close()
{
  if (m_state == State::OPEN)
    {
      if (m_output_file)
        {
          fflush (m_output_file);
          if (ferror (m_output_file))
            return Error ("error during flush");
        }

      if (m_close_file && m_output_file)
        {
          if (fclose (m_output_file) != 0)
            return Error ("error during close");

          m_output_file = nullptr;
          m_close_file = false;
        }

      m_state = State::CLOSED;
    }
  return Error::Code::NONE;
}
