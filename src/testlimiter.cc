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
using std::max;
using std::min;

class Limiter
{
  float mx = 1;
  double decay_coeff;
  uint look_ahead = 0;
  vector<float> max_buffer;
  vector<float> buffer;
public:
  Limiter (int sample_rate)
  {
    look_ahead = sample_rate * 0.005;
    assert (look_ahead >= 1);

    decay_coeff = exp (log (0.5) / (sample_rate * 0.05));
  }
  vector<float>
  process (const vector<float>& samples)
  {
    for (size_t i = 0; i < samples.size(); i++)
      {
        buffer.push_back (samples[i]);
        max_buffer.push_back (1);
      }
    for (size_t i = 0; i < buffer.size(); i++)
      {
        if (fabs (buffer[i]) > 1)
          {
            for (uint j = 0; j < look_ahead; j++)
              {
                if (int (i) - int (j) > 0)
                  {
                    double alpha = double (j) / look_ahead;
                    max_buffer[i - j] = max<float> (max_buffer[i - j], fabs (buffer[i]) * (1 - alpha) + 1 * alpha);
                  }
              }
          }
      }

    vector<float> out;
    if (buffer.size() > look_ahead)
      {
        size_t todo = buffer.size() - look_ahead;
        for (size_t i = 0; i < todo; i++)
          {
            mx = mx * decay_coeff + max_buffer[i] * (1 - decay_coeff);
            if (mx < max_buffer[i])
              mx = max_buffer[i];

            out.push_back (buffer[i] / mx);
            printf ("%f %f\n", buffer[i], out.back());
          }

        buffer.erase (buffer.begin(), buffer.begin() + todo);
        max_buffer.erase (max_buffer.begin(), max_buffer.begin() + todo);
      }
    return out;
  }
};

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
  err = out.open (in.n_channels(), in.sample_rate(), 16, in.n_frames());
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
      //out.write_frames (samples);
    }
  while (samples.size());
}
