#include <string>
#include <vector>
#include <sndfile.h>
#include <assert.h>
#include <math.h>

#include "sfinputstream.hh"
#include "stdoutwavoutputstream.hh"
#include "utils.hh"

using std::string;
using std::vector;

int
main (int argc, char **argv)
{
  SFInputStream in;
  StdoutWavOutputStream out;

  std::string filename = (argc >= 2) ? argv[1] : "-";
  Error err = in.open (filename.c_str());
  if (err)
    {
      fprintf (stderr, "teststream: open input failed: %s\n", err.message());
      return 1;
    }
  if (!out.open (in.n_channels(), in.sample_rate(), 16, in.n_frames()))
    {
      fprintf (stderr, "teststream: open output failed: %s\n", out.error_blurb());
      return 1;
    }
  vector<float> samples;
  do
    {
      in.read_frames (samples, 1024);
      out.write_frames (samples);
    }
  while (samples.size());
}
