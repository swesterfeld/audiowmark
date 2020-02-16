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
Limiter::set_block_size_ms (int ms)
{
  block_size = sample_rate * ms / 1000;
}

void
Limiter::set_ceiling (float new_ceiling)
{
  ceiling = new_ceiling;
}

vector<float>
Limiter::process (const vector<float>& samples)
{
  assert (block_size >= 1);

  const size_t n_frames = samples.size() / n_channels;
  assert (n_frames * n_channels == samples.size());    // need all channels of each frame

  buffer.insert (buffer.end(), samples.begin(), samples.end());

  /* need at least two complete blocks in buffer to produce output */
  if (buffer.size() / n_channels < 2 * block_size)
    return {};

  vector<float> out (block_size * n_channels);

  process_block (&buffer[0], &out[0]);
  buffer.erase (buffer.begin(), buffer.begin() + block_size * n_channels);

  return out;
}

void
Limiter::process_block (const float *in, float *out)
{
  float block_max = ceiling;
  float block_max2 = ceiling;
  for (size_t x = 0; x < block_size * n_channels; x++)
    {
      block_max = max (block_max, fabs (in[x]));
      block_max2 = max (block_max2, fabs (in[x + block_size * n_channels]));
    }
  float left_max = max (last_block_max, block_max);
  float right_max = max (block_max, block_max2);
  for (size_t i = 0; i < block_size; i++)
    {
      const float alpha = float (i) / block_size;
      const float scale = ceiling / (left_max * (1 - alpha) + right_max * alpha);

      for (uint c = 0; c < n_channels; c++)
        out[i * n_channels + c] = in[i * n_channels + c] * scale;
    }
  last_block_max = block_max;
}
