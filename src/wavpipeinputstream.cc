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

#include "wavpipeinputstream.hh"
#include "rawconverter.hh"

#include <assert.h>
#include <string.h>
#include <errno.h>

using std::string;
using std::vector;
using std::min;
using std::max;

WavPipeInputStream::~WavPipeInputStream()
{
  close();
}

string
header_get_4cc (unsigned char *bytes)
{
  string s;
  s += bytes[0];
  s += bytes[1];
  s += bytes[2];
  s += bytes[3];
  return s;
}

static uint16_t
header_get_u16 (unsigned char *bytes)
{
  return bytes[0] + (bytes[1] << 8);
}

static uint32_t
header_get_u32 (unsigned char *bytes)
{
  return bytes[0] + (bytes[1] << 8) + (bytes[2] << 16) + (bytes[3] << 24);
}

Error
WavPipeInputStream::read_error (const std::string& message)
{
  if (ferror (m_input_file))
    return Error (string_printf ("wav input read error: %s", strerror (errno)));
  else
    return Error (message);
}

Error
WavPipeInputStream::open (const string& filename)
{
  assert (m_state == State::NEW);

  Error err = Error::Code::NONE;
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
  RawFormat format;
  unsigned char riff_buffer[12];
  bool riff_ok = fread (riff_buffer, sizeof (riff_buffer), 1, m_input_file);
  if (!riff_ok || (header_get_4cc (riff_buffer) != "RIFF" && header_get_4cc (riff_buffer) != "RF64") ||
      header_get_4cc (riff_buffer + 8) != "WAVE")
    return read_error ("input file is not a valid wav file");

  bool in_data_chunk = false;
  bool have_fmt_chunk = false;
  do
    {
      unsigned char chunk[8];
      if (!fread (chunk, sizeof (chunk), 1, m_input_file))
        return read_error ("wav input is incomplete (no data chunk found)");

      uint32_t chunk_size = header_get_u32 (chunk + 4);
      if (header_get_4cc (chunk) == "fmt " && chunk_size >= 16 && chunk_size <= 64 * 1024 && !have_fmt_chunk)
        {
          vector<unsigned char> buffer (chunk_size);
          if (!fread (buffer.data(), buffer.size(), 1, m_input_file))
            return read_error ("wav input is incomplete (error reading fmt chunk)");

          format.set_channels (header_get_u16 (&buffer[2]));
          format.set_sample_rate (header_get_u32 (&buffer[4]));
          format.set_bit_depth (header_get_u16 (&buffer[14]));

          have_fmt_chunk = true;
        }
      else if (header_get_4cc (chunk) == "data")
        {
          in_data_chunk = true;
        }
      else // skip unknown chunk
        {
          char junk[1024];
          while (chunk_size)
            {
              uint32_t todo = min<uint32_t> (chunk_size, sizeof (junk));
              if (!fread (junk, todo, 1, m_input_file))
                return read_error ("wav input is incomplete (error skipping unknown chunk)");
              chunk_size -= todo;
            }
        }
    } while (!in_data_chunk);

  if (!have_fmt_chunk)
    return Error ("wav input is incomplete (missing fmt chunk)");

  m_raw_converter.reset (RawConverter::create (format, err));
  if (err)
    return err;

  m_format = format;
  m_state  = State::OPEN;
  return Error::Code::NONE;
}

int
WavPipeInputStream::sample_rate() const
{
  return m_format.sample_rate();
}

int
WavPipeInputStream::bit_depth() const
{
  return m_format.bit_depth();
}

size_t
WavPipeInputStream::n_frames() const
{
  return N_FRAMES_UNKNOWN;
}

int
WavPipeInputStream::n_channels() const
{
  return m_format.n_channels();
}

Error
WavPipeInputStream::read_frames (vector<float>& samples, size_t count)
{
  assert (m_state == State::OPEN);

  const size_t block_size = 1024;
  const int n_channels   = m_format.n_channels();
  const int sample_width = m_format.bit_depth() / 8;

  m_input_bytes.resize (block_size * n_channels * sample_width);
  size_t pos = 0;

  while (size_t todo = min (count, block_size))
    {
      size_t r_count = fread (m_input_bytes.data(), n_channels * sample_width, todo, m_input_file);
      if (ferror (m_input_file))
        return Error (string_printf ("error reading wav input sample data: %s", strerror (errno)));

      if (!r_count)
        break;

      samples.resize (max (samples.size(), (pos + r_count) * n_channels));

      m_raw_converter->from_raw (m_input_bytes.data(), samples.data() + pos * n_channels, r_count * n_channels);

      pos += r_count;
      count -= r_count;
    }
  samples.resize (pos * n_channels);
  return Error::Code::NONE;
}

void
WavPipeInputStream::close()
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
