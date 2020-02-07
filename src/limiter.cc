#include "limiter.hh"

#include <assert.h>
#include <math.h>
#include <stdio.h>

using std::vector;
using std::max;

Limiter::Limiter (int n_channels, int sample_rate) :
  n_channels (n_channels),
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

void
Limiter::set_ceiling (double new_ceiling)
{
  ceiling = new_ceiling;
  maximum = ceiling;
}

vector<float>
Limiter::process (const vector<float>& samples)
{
  assert (look_ahead >= 1);
  assert (release_factor > 0 && release_factor < 1);

  const size_t n_frames = samples.size() / n_channels;
  assert (n_frames * n_channels == samples.size());    // need all channels of each frame

  for (size_t i = 0; i < n_frames; i++)
    {
      for (uint c = 0; c < n_channels; c++)
        buffer.push_back (samples[i * n_channels + c]);

      max_buffer.push_back (ceiling);
    }
  for (size_t i = 0; i < n_frames; i++)
    {
      float channel_max = 0;
      for (uint c = 0; c < n_channels; c++)
        channel_max = max (channel_max, fabs (samples[i * n_channels + c]));

      if (channel_max > ceiling)
        {
          for (uint j = 0; j < look_ahead; j++)
            {
              if (int (i) - int (j) >= 0)
                {
                  double alpha = double (j) / look_ahead;
                  max_buffer[i - j] = max<float> (max_buffer[i - j], channel_max * (1 - alpha) + ceiling * alpha);
                }
            }
        }
    }

  vector<float> out;
  if (max_buffer.size() > look_ahead)
    {
      size_t todo = max_buffer.size() - look_ahead;
      for (size_t i = 0; i < todo; i++)
        {
          maximum = maximum * release_factor + max_buffer[i] * (1 - release_factor);
          if (maximum < max_buffer[i])
            maximum = max_buffer[i];

          for (uint c = 0; c < n_channels; c++)
            out.push_back (buffer[i * n_channels + c] / maximum * ceiling);
          //printf ("%f %f\n", buffer[i], out.back());
        }

      buffer.erase (buffer.begin(), buffer.begin() + todo);
      max_buffer.erase (max_buffer.begin(), max_buffer.begin() + todo);
    }
  return out;
}
