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

float
Limiter::block_max (const float *in)
{
  float maximum = ceiling;
  for (uint x = 0; x < block_size * n_channels; x++)
    maximum = max (maximum, fabs (in[x]));
  return maximum;
}

void
Limiter::process_block (const float *in, float *out)
{
  if (block_max_last < ceiling)
    block_max_last = ceiling;
  if (block_max_current < ceiling)
    block_max_current = block_max (in);
  if (block_max_next < ceiling)
    block_max_next = block_max (in + block_size * n_channels);

  const float scale_start = ceiling / max (block_max_last, block_max_current);
  const float scale_end = ceiling / max (block_max_current, block_max_next);
  const float scale_step = (scale_end - scale_start) / block_size;
  for (size_t i = 0; i < block_size; i++)
    {
      const float scale = scale_start + i * scale_step;

      // debug_scale (scale);
      for (uint c = 0; c < n_channels; c++)
        out[i * n_channels + c] = in[i * n_channels + c] * scale;
    }

  block_max_last = block_max_current;
  block_max_current = block_max_next;
  block_max_next = 0;
}

void
Limiter::debug_scale (float scale)
{
  static int debug_scale_samples = 0;

  if (debug_scale_samples % (sample_rate / 1000) == 0)
    printf ("%f %f\n", double (debug_scale_samples) / sample_rate, scale);

  debug_scale_samples++;
}

vector<float>
Limiter::flush()
{
  vector<float> out;
  vector<float> zblock (1024 * n_channels);

  size_t todo = buffer.size();
  while (todo > 0)
    {
      vector<float> block = process (zblock);
      if (block.size() > todo)
        block.resize (todo);
      out.insert (out.end(), block.begin(), block.end());
      todo -= block.size();
    }
  return out;
}
