#include "limiter.hh"

#include <assert.h>
#include <math.h>
#include <stdio.h>

using std::vector;
using std::max;

Limiter::Limiter (int sample_rate) :
  sample_rate (sample_rate)
{
}

void
Limiter::set_attack (double attack_ms)
{
  look_ahead = sample_rate / 1000.0 * attack_ms;
  look_ahead = max (look_ahead, 1u);
}

void
Limiter::set_release (double release_ms)
{
  release_factor = exp (log (0.5) / (sample_rate / 1000.0 * release_ms));
  release_factor = max (release_factor, 0.5);
}

vector<float>
Limiter::process (const vector<float>& samples)
{
  assert (look_ahead >= 1);
  assert (release_factor > 0 && release_factor < 1);

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
          maximum = maximum * release_factor + max_buffer[i] * (1 - release_factor);
          if (maximum < max_buffer[i])
            maximum = max_buffer[i];

          out.push_back (buffer[i] / maximum);
          //printf ("%f %f\n", buffer[i], out.back());
        }

      buffer.erase (buffer.begin(), buffer.begin() + todo);
      max_buffer.erase (max_buffer.begin(), max_buffer.begin() + todo);
    }
  return out;
}
