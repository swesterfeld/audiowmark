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

#include "mp3inputstream.hh"

#include <mpg123.h>
#include <assert.h>

using std::min;
using std::string;

static void
mp3_init()
{
  static bool mpg123_init_ok = false;
  if (!mpg123_init_ok)
    {
      int err = mpg123_init();
      if (err != MPG123_OK)
        {
          error ("audiowmark: init mpg123 lib failed\n");
          exit (1);
        }
      mpg123_init_ok = true;
    }
}

MP3InputStream::~MP3InputStream()
{
  close();
}

Error
MP3InputStream::open (const string& filename)
{
  int err = 0;

  mp3_init();

  m_handle = mpg123_new (nullptr, &err);
  if (err != MPG123_OK)
    return Error ("mpg123_new failed");

  err = mpg123_param (m_handle, MPG123_ADD_FLAGS, MPG123_QUIET, 0);
  if (err != MPG123_OK)
    return Error ("setting quiet mode failed");

  // allow arbitary amount of data for resync */
  err = mpg123_param (m_handle, MPG123_RESYNC_LIMIT, -1, 0);
  if (err != MPG123_OK)
    return Error ("setting resync limit parameter failed");

  // force floating point output
  {
    const long *rates;
    size_t      rate_count;

    mpg123_format_none (m_handle);
    mpg123_rates (&rates, &rate_count);

    for (size_t i = 0; i < rate_count; i++)
      {
        err = mpg123_format (m_handle, rates[i], MPG123_MONO|MPG123_STEREO, MPG123_ENC_FLOAT_32);
        if (err != MPG123_OK)
          return Error (mpg123_strerror (m_handle));
      }
  }

  err = mpg123_open (m_handle, filename.c_str());
  if (err != MPG123_OK)
    return Error (mpg123_strerror (m_handle));

  m_need_close = true;

  /* scan headers to get best possible length estimate */
  err = mpg123_scan (m_handle);
  if (err != MPG123_OK)
    return Error (mpg123_strerror (m_handle));

  long rate;
  int channels;
  int encoding;

  err = mpg123_getformat (m_handle, &rate, &channels, &encoding);
  if (err != MPG123_OK)
    return Error (mpg123_strerror (m_handle));

  /* ensure that the format will not change */
  mpg123_format_none (m_handle);
  mpg123_format (m_handle, rate, channels, encoding);

  m_n_values = mpg123_length (m_handle) * channels;
  m_n_channels = channels;
  m_sample_rate = rate;
  m_frames_left = m_n_values / m_n_channels;

  return Error::Code::NONE;
}

Error
MP3InputStream::read_frames (std::vector<float>& samples, size_t count)
{
  while (!m_eof && m_read_buffer.size() < count * m_n_channels)
    {
      size_t buffer_bytes = mpg123_outblock (m_handle);
      assert (buffer_bytes % sizeof (float) == 0);

      std::vector<float> buffer (buffer_bytes / sizeof (float));

      size_t done;
      int err = mpg123_read (m_handle, reinterpret_cast<unsigned char *> (&buffer[0]), buffer_bytes, &done);
      if (err == MPG123_OK)
        {
          const size_t n_values = done / sizeof (float);
          m_read_buffer.insert (m_read_buffer.end(), buffer.begin(), buffer.begin() + n_values);
        }
      else if (err == MPG123_DONE)
        {
          m_eof = true;
        }
      else if (err == MPG123_NEED_MORE)
        {
          // some mp3s have this error before reaching eof -> harmless
          m_eof = true;
        }
      else
        {
          return Error (mpg123_strerror (m_handle));
        }
    }
  /* pad zero samples at end if necessary to match the number of frames we promised to deliver */
  if (m_eof && m_read_buffer.size() < m_frames_left * m_n_channels)
    m_read_buffer.resize (m_frames_left * m_n_channels);

  /* never read past the promised number of frames */
  if (count > m_frames_left)
    count = m_frames_left;

  const auto begin = m_read_buffer.begin();
  const auto end   = begin + min (count * m_n_channels, m_read_buffer.size());
  samples.assign (begin, end);
  m_read_buffer.erase (begin, end);
  m_frames_left -= count;
  return Error::Code::NONE;
}

void
MP3InputStream::close()
{
  if (m_state == State::OPEN)
    {
      if (m_handle && m_need_close)
        mpg123_close (m_handle);

      if (m_handle)
        {
          mpg123_delete (m_handle);
          m_handle = nullptr;
        }
      m_state = State::CLOSED;
    }
}

int
MP3InputStream::bit_depth() const
{
  return 24; /* mp3 decoder is running on floats */
}

int
MP3InputStream::sample_rate() const
{
  return m_sample_rate;
}

int
MP3InputStream::n_channels() const
{
  return m_n_channels;
}

size_t
MP3InputStream::n_frames() const
{
  return m_n_values / m_n_channels;
}

/* there is no really simple way of detecting if something is an mp3
 *
 * so we try to decode a few frames; if that works without error the
 * file is probably a valid mp3
 */
bool
MP3InputStream::detect (const string& filename)
{
  struct ScopedMHandle
  {
    mpg123_handle *mh         = nullptr;
    bool           need_close = false;

    ~ScopedMHandle()
    {
      if (mh && need_close)
        mpg123_close (mh);

      if (mh)
        mpg123_delete (mh);
    }
  };

  int err = 0;

  mp3_init();

  mpg123_handle *mh = mpg123_new (NULL, &err);
  if (err != MPG123_OK)
    return false;

  auto smh = ScopedMHandle { mh }; // cleanup on return

  err = mpg123_param (mh, MPG123_ADD_FLAGS, MPG123_QUIET, 0);
  if (err != MPG123_OK)
    return false;

  err = mpg123_open (mh, filename.c_str());
  if (err != MPG123_OK)
    return false;

  smh.need_close = true;

  long rate;
  int channels;
  int encoding;
  err = mpg123_getformat (mh, &rate, &channels, &encoding);
  if (err != MPG123_OK)
    return false;

  size_t buffer_bytes = mpg123_outblock (mh);
  unsigned char buffer[buffer_bytes];

  for (size_t i = 0; i < 30; i++)
    {
      size_t done;
      err = mpg123_read (mh, buffer, buffer_bytes, &done);
      if (err == MPG123_DONE)
        {
          return true;
        }
      else if (err != MPG123_OK)
        {
          return false;
        }
    }
  return true;
}
