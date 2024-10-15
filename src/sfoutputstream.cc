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

#include "sfoutputstream.hh"
#include "rawconverter.hh"
#include "utils.hh"
#include "wmcommon.hh"

#include <math.h>
#include <assert.h>

using std::string;
using std::vector;

SFOutputStream::~SFOutputStream()
{
  close();
}

Error
SFOutputStream::open (const string& filename, int n_channels, int sample_rate, int bit_depth, OutFormat out_format)
{
  return open ([&] (SF_INFO *sfinfo) {
    return sf_open (filename.c_str(), SFM_WRITE, sfinfo);
  }, n_channels, sample_rate, bit_depth, out_format);
}

// supported WAV subtypes
static constexpr uint64_t WAV_SUBFORMATS =
  SF_FORMAT_PCM_U8    |
  SF_FORMAT_PCM_16    |
  SF_FORMAT_PCM_24    |
  SF_FORMAT_PCM_32    |
  SF_FORMAT_FLOAT     |
  SF_FORMAT_DOUBLE    |
  SF_FORMAT_ULAW      |
  SF_FORMAT_ALAW      |
  SF_FORMAT_IMA_ADPCM |
  SF_FORMAT_MS_ADPCM  |
  0;


Error
SFOutputStream::open (std::function<SNDFILE* (SF_INFO *)> open_func, int n_channels, int sample_rate, int bit_depth, OutFormat out_format)
{
  assert (m_state == State::NEW);

  m_sample_rate = sample_rate;
  m_n_channels  = n_channels;

  SF_INFO sfinfo = {0,};
  sfinfo.samplerate = sample_rate;
  sfinfo.channels   = n_channels;

  if (bit_depth > 16)
    {
      sfinfo.format = SF_FORMAT_PCM_24;
      m_bit_depth   = 24;
    }
  else
    {
      sfinfo.format = SF_FORMAT_PCM_16;
      m_bit_depth   = 16;
    }
   switch (out_format)
     {
       case OutFormat::WAV:
         sfinfo.format |= SF_FORMAT_WAV;
         if (Params::output_format == Format::AUTO && (WAV_SUBFORMATS & Params::input_sndfile_flags)) {
           sfinfo.format = SF_FORMAT_WAV | (WAV_SUBFORMATS & Params::input_sndfile_flags);      /* preserve input format */
           m_write_float_data = (sfinfo.format & SF_FORMAT_DOUBLE) == SF_FORMAT_DOUBLE ||
                                (sfinfo.format & SF_FORMAT_FLOAT) == SF_FORMAT_FLOAT;
         }
         break;
       case OutFormat::RF64: sfinfo.format |= SF_FORMAT_RF64;
                             break;
       case OutFormat::FLAC: sfinfo.format |= SF_FORMAT_FLAC;
                             break;
       default:              assert (false);
     }

  m_sndfile = open_func (&sfinfo);
  int error = sf_error (m_sndfile);
  if (error)
    {
      string msg = sf_strerror (m_sndfile);
      if (m_sndfile)
        sf_close (m_sndfile);

      return Error (msg);
    }
  m_state       = State::OPEN;
  return Error::Code::NONE;
}

Error
SFOutputStream::close()
{
  if (m_state == State::OPEN)
    {
      assert (m_sndfile);
      if (sf_close (m_sndfile))
        return Error ("sf_close returned an error");

      m_state = State::CLOSED;
    }
  return Error::Code::NONE;
}

Error
SFOutputStream::write_frames (const vector<float>& samples)
{
  const sf_count_t frames = samples.size() / m_n_channels;
  sf_count_t count;
  if (m_write_float_data)
    count = sf_writef_float (m_sndfile, samples.data(), frames);
  else {
    vector<int> isamples (samples.size());
    for (size_t i = 0; i < samples.size(); i++)
      isamples[i] = float_to_int_clip<32> (samples[i]);
    count = sf_writef_int (m_sndfile, isamples.data(), frames);
  }

  if (sf_error (m_sndfile))
    return Error (sf_strerror (m_sndfile));

  if (count != frames)
    return Error ("writing sample data failed: short write");

  return Error::Code::NONE;
}

int
SFOutputStream::bit_depth() const
{
  return m_bit_depth;
}

int
SFOutputStream::sample_rate() const
{
  return m_sample_rate;
}

int
SFOutputStream::n_channels() const
{
  return m_n_channels;
}

Error
SFOutputStream::open (vector<unsigned char> *data, int n_channels, int sample_rate, int bit_depth, OutFormat out_format)
{
  m_virtual_data.mem = data;

  return open ([&] (SF_INFO *sfinfo) {
    return sf_open_virtual (&m_virtual_data.io, SFM_WRITE, sfinfo, &m_virtual_data);
  }, n_channels, sample_rate, bit_depth, out_format);
}
