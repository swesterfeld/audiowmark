#include "audiostream.hh"
#include "wmcommon.hh"
#include "sfinputstream.hh"
#include "sfoutputstream.hh"
#include "mp3.hh"
#include "mp3inputstream.hh"
#include "rawconverter.hh"
#include "rawoutputstream.hh"
#include "stdoutwavoutputstream.hh"

using std::string;

AudioStream::~AudioStream()
{
}

std::unique_ptr<AudioInputStream>
AudioInputStream::create (const string& filename, Error& err)
{
  std::unique_ptr<AudioInputStream> in_stream;

  if (Params::input_format == Format::AUTO)
    {
      SFInputStream *sistream = new SFInputStream();
      in_stream.reset (sistream);
      err = sistream->open (filename);
      if (err && MP3InputStream::detect (filename))
        {
          MP3InputStream *mistream = new MP3InputStream();
          in_stream.reset (mistream);

          err = mistream->open (filename);
          if (err)
            return nullptr;
        }
      else if (err)
        return nullptr;
    }
  else
    {
      RawInputStream *ristream = new RawInputStream();
      in_stream.reset (ristream);

      err = ristream->open (filename, Params::raw_input_format);
      if (err)
        return nullptr;
    }
  return in_stream;
}

std::unique_ptr<AudioOutputStream>
AudioOutputStream::create (const string& filename, int n_channels, int sample_rate, int bit_depth, size_t n_frames, Error& err)
{
  std::unique_ptr<AudioOutputStream> out_stream;

  if (Params::output_format == Format::RAW)
    {
      RawOutputStream *rostream = new RawOutputStream();
      out_stream.reset (rostream);
      err = rostream->open (filename, Params::raw_output_format);
      if (err)
        return nullptr;
    }
  else if (filename == "-")
    {
      StdoutWavOutputStream *swstream = new StdoutWavOutputStream();
      out_stream.reset (swstream);
      err = swstream->open (n_channels, sample_rate, bit_depth, n_frames);
      if (err)
        return nullptr;
    }
  else
    {
      SFOutputStream *sfostream = new SFOutputStream();
      out_stream.reset (sfostream);
      err = sfostream->open (filename, n_channels, sample_rate, bit_depth, n_frames);
      if (err)
        return nullptr;
    }
  return out_stream;
}
