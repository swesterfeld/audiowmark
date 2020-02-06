#include <string>
#include <vector>
#include <sndfile.h>
#include <assert.h>
#include <math.h>

#include "sfinputstream.hh"
#include "sfoutputstream.hh"
#include "utils.hh"
#include "limiter.hh"

using std::string;
using std::vector;
using std::max;
using std::min;


int
main (int argc, char **argv)
{
  SFInputStream in;
  SFOutputStream out;

  Error err = in.open (argv[1]);
  if (err)
    {
      fprintf (stderr, "teststream: open input failed: %s\n", err.message());
      return 1;
    }
  err = out.open (argv[2], in.n_channels(), in.sample_rate(), 16, in.n_frames());
  if (err)
    {
      fprintf (stderr, "teststream: open output failed: %s\n", err.message());
      return 1;
    }
  Limiter limiter (in.sample_rate());
  vector<float> samples;
  do
    {
      in.read_frames (samples, 1024);
      for (auto& s: samples)
        s *= 1.1;
      samples = limiter.process (samples);
      out.write_frames (samples);
    }
  while (samples.size());
}
