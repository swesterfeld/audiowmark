#ifndef AUDIOWMARK_LIMITER_HH
#define AUDIOWMARK_LIMITER_HH

#include <vector>
#include <sys/types.h>

class Limiter
{
  float ceiling         = 1;
  float last_block_max  = 0;
  uint  block_size      = 0;
  uint  n_channels      = 0;
  uint  sample_rate     = 0;

  std::vector<float> buffer;
  void process_block (const float *in, float *out);
  void debug_scale (float scale);
public:
  Limiter (int n_channels, int sample_rate);

  void set_block_size_ms (int value_ms);
  void set_ceiling (float ceiling);

  std::vector<float> process (const std::vector<float>& samples);
};

#endif /* AUDIOWMARK_LIMITER_HH */
