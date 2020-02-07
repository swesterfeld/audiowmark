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
      fprintf (stderr, "testlimiter: open input failed: %s\n", err.message());
      return 1;
    }
  err = out.open (argv[2], in.n_channels(), in.sample_rate(), 16, in.n_frames());
  if (err)
    {
      fprintf (stderr, "testlimiter: open output failed: %s\n", err.message());
      return 1;
    }
  Limiter limiter (in.n_channels(), in.sample_rate());
  limiter.set_attack (5);
  limiter.set_release (50);
  limiter.set_ceiling (0.9);
  vector<float> in_samples;
  do
    {
      in.read_frames (in_samples, 1024);
      for (auto& s: in_samples)
        s *= 1.1;

      vector<float> out_samples = limiter.process (in_samples);
      out.write_frames (out_samples);
    }
  while (in_samples.size());
}
