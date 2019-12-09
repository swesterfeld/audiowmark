#include "wavdata.hh"
#include "mp3.hh"
#include "utils.hh"
#include "audiostream.hh"
#include "sfinputstream.hh"
#include "sfoutputstream.hh"
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
  std::unique_ptr<AudioOutputStream> out_stream; // FIXME: virtual constructor

  SFOutputStream *sostream = new SFOutputStream();
  out_stream.reset (sostream);
  Error err = sostream->open (filename, m_n_channels, m_sample_rate, m_bit_depth, m_samples.size() / m_n_channels);
  if (err)
    {
      m_error_blurb = err.message();
      return false;
    }

  err = sostream->write_frames (m_samples);
  if (err)
    {
      m_error_blurb = err.message();
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
