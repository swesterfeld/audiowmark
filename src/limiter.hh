#ifndef AUDIOWMARK_LIMITER_HH
#define AUDIOWMARK_LIMITER_HH

#include <vector>
#include <sys/types.h>

class Limiter
{
  float ceiling           = 1;
  float block_max_last    = 0;
  float block_max_current = 0;
  float block_max_next    = 0;
  uint  block_size        = 0;
  uint  n_channels        = 0;
  uint  sample_rate       = 0;
  size_t buffered_frames  = 0;

  std::vector<float> buffer;
  void process_block (const float *in, float *out);
  float block_max (const float *in);
  void debug_scale (float scale);
public:
  Limiter (int n_channels, int sample_rate);

  void set_block_size_ms (int value_ms);
  void set_ceiling (float ceiling);

  std::vector<float> process (const std::vector<float>& samples);
  std::vector<float> flush();
};

#endif /* AUDIOWMARK_LIMITER_HH */
