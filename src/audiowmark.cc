#include <string.h>
#include <string>

#include "wavdata.hh"

using std::string;

int
add_watermark (const string& infile, const string& outfile, const string& bits)
{
  printf ("loading %s\n", infile.c_str());

  WavData wav_data;
  if (!wav_data.load (infile))
    {
      fprintf (stderr, "audiowmark: error loading %s: %s\n", infile.c_str(), wav_data.error_blurb());
      return 1;
    }
  printf ("channels: %d, samples: %zd, mix_freq: %f\n", wav_data.n_channels(), wav_data.n_values(), wav_data.mix_freq());

  // magic in here...

  if (!wav_data.save (outfile))
    {
      fprintf (stderr, "audiowmark: error saving %s: %s\n", outfile.c_str(), wav_data.error_blurb());
      return 1;
    }
}

int
main (int argc, char **argv)
{
  if (strcmp (argv[1], "add") == 0 && argc == 5)
    {
      return add_watermark (argv[2], argv[3], argv[4]);
    }
  else
    {
      fprintf (stderr, "audiowmark: error parsing commandline args\n");
      return 1;
    }
}
