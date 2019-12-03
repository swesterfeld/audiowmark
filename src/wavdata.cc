#include "wavdata.hh"
#include "mp3.hh"
#include "utils.hh"
#include "audiostream.hh"
#include "sfinputstream.hh"
#include "mp3inputstream.hh"

#include <memory>
#include <math.h>
#include <sndfile.h>

using std::string;
using std::vector;

WavData::WavData()
{
}

WavData::WavData (const vector<float>& samples, int n_channels, int sample_rate, int bit_depth)
{
  m_samples     = samples;
  m_n_channels  = n_channels;
  m_sample_rate = sample_rate;
  m_bit_depth   = bit_depth;
}

Error
WavData::load (const string& filename)
{
  std::unique_ptr<AudioInputStream> in_stream; // FIXME: virtual constructor

  SFInputStream *sistream = new SFInputStream();
  in_stream.reset (sistream);
  Error err = sistream->open (filename);
  if (err && mp3_detect (filename))
    {
      MP3InputStream *mistream = new MP3InputStream();
      in_stream.reset (mistream);

      err = mistream->open (filename);
      if (err)
        return err;
    }
  else if (err)
    return err;

  vector<float> m_buffer;
  while (true)
    {
      err = in_stream->read_frames (m_buffer, 1024);
      if (err)
        return err;

      if (!m_buffer.size())
        {
          /* reached eof */
          break;
        }
      m_samples.insert (m_samples.end(), m_buffer.begin(), m_buffer.end());
    }
  m_sample_rate = in_stream->sample_rate();
  m_n_channels  = in_stream->n_channels();
  m_bit_depth   = in_stream->bit_depth();

  return Error::Code::NONE;
}

bool
WavData::save (const string& filename)
{
  SF_INFO sfinfo = {0,};

  sfinfo.samplerate = m_sample_rate;
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

int
WavData::sample_rate() const
{
  return m_sample_rate;
}

int
WavData::bit_depth() const
{
  return m_bit_depth;
}

void
WavData::set_samples (const vector<float>& samples)
{
  m_samples = samples;
}

const char *
WavData::error_blurb() const
{
  return m_error_blurb.c_str();
}
