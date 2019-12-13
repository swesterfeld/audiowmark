#include "audiostream.hh"
#include "wmcommon.hh"
#include "sfinputstream.hh"
#include "mp3.hh"
#include "mp3inputstream.hh"
#include "rawconverter.hh"

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
      if (err && mp3_detect (filename))
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
