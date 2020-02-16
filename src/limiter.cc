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
  const uint buffered_blocks = buffer.size() / n_channels / block_size;
  if (buffered_blocks < 2)
    return {};

  const uint blocks_todo = buffered_blocks - 1;

  vector<float> out (blocks_todo * block_size * n_channels);
  for (uint b = 0; b < blocks_todo; b++)
    process_block (&buffer[b * block_size * n_channels], &out[b * block_size * n_channels]);

  buffer.erase (buffer.begin(), buffer.begin() + blocks_todo * block_size * n_channels);

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
  const float scale_start = ceiling / max (last_block_max, block_max);
  const float scale_end = ceiling / max (block_max, block_max2);
  for (size_t i = 0; i < block_size; i++)
    {
      const float alpha = float (i) / block_size;
      const float scale = scale_start * (1 - alpha) + scale_end * alpha;

      // debug_scale (scale);
      for (uint c = 0; c < n_channels; c++)
        out[i * n_channels + c] = in[i * n_channels + c] * scale;
    }
  last_block_max = block_max;
}

void
Limiter::debug_scale (float scale)
{
  static int debug_scale_samples = 0;

  if (debug_scale_samples % (sample_rate / 1000) == 0)
    printf ("%f %f\n", double (debug_scale_samples) / sample_rate, scale);

  debug_scale_samples++;
}
