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

#include "sfinputstream.hh"

#include <assert.h>
#include <string.h>
#include <unistd.h>

using std::string;
using std::vector;

SFInputStream::~SFInputStream()
{
  close();
}

Error
SFInputStream::open (const string& filename)
{
  return open ([&] (SF_INFO *sfinfo) {
    if (filename == "-")
      {
        m_is_stdin = true;
        return sf_open_fd (STDIN_FILENO, SFM_READ, sfinfo, /* close fd */ SF_FALSE);
      }
    else
      {
        return sf_open (filename.c_str(), SFM_READ, sfinfo);
      }
  });
}


Error
SFInputStream::open (std::function<SNDFILE* (SF_INFO *)> open_func)
{
  assert (m_state == State::NEW);

  SF_INFO sfinfo = { 0, };

  m_sndfile = open_func (&sfinfo);

  int error = sf_error (m_sndfile);
  if (error)
    {
      Error err (sf_strerror (m_sndfile));
      if (m_sndfile)
        {
          m_sndfile = nullptr;
          sf_close (m_sndfile);
        }
      return err;
    }

  m_n_channels  = sfinfo.channels;
  m_n_frames    = (sfinfo.frames == SF_COUNT_MAX) ? N_FRAMES_UNKNOWN : sfinfo.frames;
  m_sample_rate = sfinfo.samplerate;

  switch (sfinfo.format & SF_FORMAT_SUBMASK)
    {
      case SF_FORMAT_PCM_U8:
      case SF_FORMAT_PCM_S8:
          m_bit_depth = 8;
          break;

      case SF_FORMAT_PCM_16:
          m_bit_depth = 16;
          break;

      case SF_FORMAT_PCM_24:
          m_bit_depth = 24;
          break;

      case SF_FORMAT_PCM_32:
          m_bit_depth = 32;
          break;

      case SF_FORMAT_FLOAT:
          m_bit_depth = 32;
          m_read_float_data = true;
          break;

      case SF_FORMAT_DOUBLE:
          m_bit_depth = 64;
          m_read_float_data = true;
          break;

      default:
          m_bit_depth = 32; /* unknown */
    }

  m_state       = State::OPEN;
  return Error::Code::NONE;
}

int
SFInputStream::sample_rate() const
{
  return m_sample_rate;
}

int
SFInputStream::bit_depth() const
{
  return m_bit_depth;
}

Error
SFInputStream::read_frames (vector<float>& samples, size_t count)
{
  assert (m_state == State::OPEN);

  if (m_read_float_data) /* float or double input */
    {
      samples.resize (count * m_n_channels);

      sf_count_t r_count = sf_readf_float (m_sndfile, &samples[0], count);

      if (sf_error (m_sndfile))
        return Error (sf_strerror (m_sndfile));

      samples.resize (r_count * m_n_channels);
    }
  else /* integer input */
    {
      vector<int> isamples (count * m_n_channels);

      sf_count_t r_count = sf_readf_int (m_sndfile, &isamples[0], count);

      if (sf_error (m_sndfile))
        return Error (sf_strerror (m_sndfile));

      /* reading a wav file and saving it again with the libsndfile float API will
       * change some values due to normalization issues:
       *   http://www.mega-nerd.com/libsndfile/FAQ.html#Q010
       *
       * to avoid the problem, we use the int API and do the conversion beween int
       * and float manually - the important part is that the normalization factors
       * used during read and write are identical
       */
      samples.resize (r_count * m_n_channels);
      const double norm = 1.0 / 0x80000000LL;
      for (size_t i = 0; i < samples.size(); i++)
        samples[i] = isamples[i] * norm;
    }

  return Error::Code::NONE;
}

void
SFInputStream::close()
{
  if (m_state == State::OPEN)
    {
      assert (m_sndfile);
      sf_close (m_sndfile);

      m_sndfile = nullptr;
      m_state = State::CLOSED;

      if (m_is_stdin)
        {
          /* WAV files can contain additional RIFF chunks after the end of the 'data' chunk (issue #19).
           *  -> skip the rest of stdin to avoid SIGPIPE for the process writing to the pipe
           */
          ssize_t count;

          do
            {
              char junk[16 * 1024];
              count = read (STDIN_FILENO, junk, sizeof (junk)) ;
            }
          while (count > 0 || (count == -1 && errno == EINTR));
        }
    }
}

static sf_count_t
virtual_get_len (void *data)
{
  SFVirtualData *vdata = static_cast<SFVirtualData *> (data);

  return vdata->mem->size();
}

static sf_count_t
virtual_seek (sf_count_t offset, int whence, void *data)
{
  SFVirtualData *vdata = static_cast<SFVirtualData *> (data);

  if (whence == SEEK_CUR)
    {
      vdata->offset = vdata->offset + offset;
    }
  else if (whence == SEEK_SET)
    {
      vdata->offset = offset;
    }
  else if (whence == SEEK_END)
    {
      vdata->offset = vdata->mem->size() + offset;
    }

  /* can't seek beyond eof */
  vdata->offset = bound<sf_count_t> (0, vdata->offset, vdata->mem->size());
  return vdata->offset;
}

static sf_count_t
virtual_read (void *ptr, sf_count_t count, void *data)
{
  SFVirtualData *vdata = static_cast<SFVirtualData *> (data);

  int rcount = 0;
  if (size_t (vdata->offset + count) <= vdata->mem->size())
    {
      /* fast case: read can be fully satisfied with the data we have */
      memcpy (ptr, &(*vdata->mem)[vdata->offset], count);
      rcount = count;
    }
  else
    {
      unsigned char *uptr = static_cast<unsigned char *> (ptr);
      for (sf_count_t i = 0; i < count; i++)
        {
          size_t rpos = i + vdata->offset;
          if (rpos < vdata->mem->size())
            {
              uptr[i] = (*vdata->mem)[rpos];
              rcount++;
            }
        }
    }
  vdata->offset += rcount;
  return rcount;
}

static sf_count_t
virtual_write (const void *ptr, sf_count_t count, void *data)
{
  SFVirtualData *vdata = static_cast<SFVirtualData *> (data);

  const unsigned char *uptr = static_cast<const unsigned char *> (ptr);
  for (sf_count_t i = 0; i < count; i++)
    {
      unsigned char ch = uptr[i];

      size_t wpos = i + vdata->offset;
      if (wpos >= vdata->mem->size())
        vdata->mem->resize (wpos + 1);
      (*vdata->mem)[wpos] = ch;
    }
  vdata->offset += count;
  return count;
}

static sf_count_t
virtual_tell (void *data)
{
  SFVirtualData *vdata = static_cast<SFVirtualData *> (data);
  return vdata->offset;
}

SFVirtualData::SFVirtualData() :
  io {
    virtual_get_len,
    virtual_seek,
    virtual_read,
    virtual_write,
    virtual_tell
  }
{
}

Error
SFInputStream::open (const vector<unsigned char> *data)
{
  m_virtual_data.mem = const_cast<vector<unsigned char> *> (data);
  return open ([&] (SF_INFO *sfinfo) {
    return sf_open_virtual (&m_virtual_data.io, SFM_READ, sfinfo, &m_virtual_data);
  });
}
