#include "rawinputstream.hh"

#include <assert.h>

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
#if 0
  SF_INFO sfinfo = { 0, };

  m_sndfile = sf_open (filename.c_str(), SFM_READ, &sfinfo);

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
#endif

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
  size_t r_count = fread (input_bytes.data(), n_channels * sample_width, count, stdin);

  unsigned char *ptr = reinterpret_cast<unsigned char *> (input_bytes.data());

  samples.resize (r_count * n_channels);
  const double norm = 1.0 / 0x80000000LL;
  for (size_t i = 0; i < samples.size(); i++)
    {
      int s32 = (ptr[1] << 24) + (ptr[0] << 16);
      samples[i] = s32 * norm;
      ptr += 2;
    }

  return Error::Code::NONE;
}

void
RawInputStream::close()
{
  if (m_state == State::OPEN)
    {
#if 0
      assert (m_sndfile);
      sf_close (m_sndfile);
#endif

      m_state = State::CLOSED;
    }
}
