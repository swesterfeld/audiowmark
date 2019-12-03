#include "mp3inputstream.hh"
#include "mp3.hh"

#include <mpg123.h>
#include <assert.h>

using std::min;

MP3InputStream::~MP3InputStream()
{
  close();
}

Error
MP3InputStream::open (const std::string& filename)
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
  mpg123_scan (m_handle);

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
  m_frames_left -= samples.size() / m_n_channels;
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
