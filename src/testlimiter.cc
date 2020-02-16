#include <string>
#include <vector>
#include <sndfile.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#include "sfinputstream.hh"
#include "sfoutputstream.hh"
#include "utils.hh"
#include "limiter.hh"

using std::string;
using std::vector;
using std::max;
using std::min;

static double
gettime()
{
  timeval tv;
  gettimeofday (&tv, 0);

  return tv.tv_sec + tv.tv_usec / 1000000.0;
}

int
perf()
{
  Limiter limiter (2, 44100);

  limiter.set_block_size_ms (1000);

  vector<float> samples (2 * 1024);

  int n_frames = 0;
  double start = gettime();
  for (int i = 0; i < 100000; i++)
    {
      n_frames += samples.size() / 2;
      vector<float> out_samples = limiter.process (samples);
    }
  double end = gettime();
  printf ("%f ns/frame\n", (end - start) * 1000 * 1000 * 1000 / n_frames);
  return 0;
}

int
impulses()
{
  Limiter limiter (2, 44100);
  limiter.set_block_size_ms (3);
  limiter.set_ceiling (0.9);

  vector<float> in_all, out_all;
  int pos = 0;
  for (int block = 0; block < 10; block++)
    {
      vector<float> in_samples;
      for (int i = 0; i < 1024; i++)
        {
          double d = (pos++ % 441) == 440 ? 1.0 : 0.5;
          in_samples.push_back (d);
          in_samples.push_back (d); /* stereo */
        }
      vector<float> out_samples = limiter.process (in_samples);

      in_all.insert (in_all.end(), in_samples.begin(), in_samples.end());
      out_all.insert (out_all.end(), out_samples.begin(), out_samples.end());
    }
  for (size_t i = 0; i < min (in_all.size(), out_all.size()); i += 2)
    {
      assert (out_all[i] == out_all[i + 1]); /* stereo */
      printf ("%f %f\n", in_all[i], out_all[i]);
    }
  return 0;
}

int
main (int argc, char **argv)
{
  if (argc == 2 && strcmp (argv[1], "perf") == 0)
    return perf();
  if (argc == 2 && strcmp (argv[1], "impulses") == 0)
    return impulses();

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
  limiter.set_block_size_ms (1000);
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
