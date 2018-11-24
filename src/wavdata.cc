#include "wavdata.hh"

#include <math.h>
#include <sndfile.h>

using std::string;
using std::vector;

template<typename T>
inline const T&
bound (const T& min_value, const T& value, const T& max_value)
{
  return std::min (std::max (value, min_value), max_value);
}


bool
WavData::load (const string& filename)
{
  SF_INFO sfinfo = { 0, };

  SNDFILE *sndfile = sf_open (filename.c_str(), SFM_READ, &sfinfo);

  int error = sf_error (sndfile);
  if (error)
    {
      m_error_blurb = sf_strerror (sndfile);
      if (sndfile)
        sf_close (sndfile);

      return false;
    }

  vector<int> isamples (sfinfo.frames * sfinfo.channels);
  sf_count_t count = sf_readf_int (sndfile, &isamples[0], sfinfo.frames);

  error = sf_error (sndfile);
  if (error)
    {
      m_error_blurb = sf_strerror (sndfile);
      sf_close (sndfile);

      return false;
    }

  if (count != sfinfo.frames)
    {
      m_error_blurb = "reading sample data failed: short read";
      sf_close (sndfile);

      return false;
    }

  m_samples.resize (sfinfo.frames * sfinfo.channels);

  /* reading a wav file and saving it again with the libsndfile float API will
   * change some values due to normalization issues:
   *   http://www.mega-nerd.com/libsndfile/FAQ.html#Q010
   *
   * to avoid the problem, we use the int API and do the conversion beween int
   * and float manually - the important part is that the normalization factors
   * used during read and write are identical
   */
  const double norm = 1.0 / 0x80000000LL;
  for (size_t i = 0; i < m_samples.size(); i++)
    m_samples[i] = isamples[i] * norm;

  m_mix_freq    = sfinfo.samplerate;
  m_n_channels  = sfinfo.channels;

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

      case SF_FORMAT_FLOAT:
      case SF_FORMAT_PCM_32:
          m_bit_depth = 32;
          break;

      case SF_FORMAT_DOUBLE:
          m_bit_depth = 64;
          break;

      default:
          m_bit_depth = 32; /* unknown */
    }

  error = sf_close (sndfile);
  if (error)
    {
      m_error_blurb = sf_error_number (error);
      return false;
    }
  return true;
}

bool
WavData::save (const string& filename)
{
  SF_INFO sfinfo = {0,};

  sfinfo.samplerate = lrint (m_mix_freq);
  sfinfo.channels   = m_n_channels;

  if (m_bit_depth > 16)
    sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_24;
  else
    sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;

  SNDFILE *sndfile = sf_open (filename.c_str(), SFM_WRITE, &sfinfo);
  int error = sf_error (sndfile);
  if (error)
    {
      m_error_blurb = sf_strerror (sndfile);
      if (sndfile)
        sf_close (sndfile);

      return false;
    }

  vector<int> isamples (m_samples.size());
  for (size_t i = 0; i < m_samples.size(); i++)
    {
      const double norm      =  0x80000000LL;
      const double min_value = -0x80000000LL;
      const double max_value =  0x7FFFFFFF;

      isamples[i] = lrint (bound<double> (min_value, m_samples[i] * norm, max_value));
    }

  sf_count_t frames = m_samples.size() / m_n_channels;
  sf_count_t count = sf_writef_int (sndfile, &isamples[0], frames);

  error = sf_error (sndfile);
  if (error)
    {
      m_error_blurb = sf_strerror (sndfile);
      sf_close (sndfile);

      return false;
    }

  if (count != frames)
    {
      m_error_blurb = "writing sample data failed: short write";
      sf_close (sndfile);

      return false;
    }

  error = sf_close (sndfile);
  if (error)
    {
      m_error_blurb = sf_error_number (error);
      return false;
    }
  return true;
}

float
WavData::mix_freq() const
{
  return m_mix_freq;
}

int
WavData::n_channels() const
{
  return m_n_channels;
}

int
WavData::bit_depth() const
{
  return m_bit_depth;
}

const vector<float>&
WavData::samples() const
{
  return m_samples;
}

size_t
WavData::n_values() const
{
  return m_samples.size();
}

const char *
WavData::error_blurb() const
{
  return m_error_blurb.c_str();
}
