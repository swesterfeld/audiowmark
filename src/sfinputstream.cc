#include "sfinputstream.hh"

#include <assert.h>

using std::string;
using std::vector;

SFInputStream::~SFInputStream()
{
  close();
}

bool
SFInputStream::open (const string& filename)
{
  assert (m_state == State::NEW);

  SF_INFO sfinfo = { 0, };

  m_sndfile = sf_open (filename.c_str(), SFM_READ, &sfinfo);

  int error = sf_error (m_sndfile);
  if (error)
    {
      m_error_blurb = sf_strerror (m_sndfile);
      if (m_sndfile)
        {
          m_sndfile = nullptr;
          sf_close (m_sndfile);
        }
      return false;
    }

  m_n_channels  = sfinfo.channels;
  m_n_values    = sfinfo.frames * sfinfo.channels;
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

  m_state       = State::OPEN;
  return true;
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

  return Error::Code::NONE;
}

void
SFInputStream::close()
{
  if (m_state == State::OPEN)
    {
      assert (m_sndfile);
      sf_close (m_sndfile);

      m_state = State::CLOSED;
    }
}
