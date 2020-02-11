#ifndef AUDIOWMARK_LIMITER_HH
#define AUDIOWMARK_LIMITER_HH

#include <vector>
#include <sys/types.h>

class Limiter
{
  float ceiling         = 1;
  float maximum         = 1;
  float release_factor  = 0;
  uint look_ahead       = 0;
  uint n_channels       = 0;
  uint sample_rate      = 0;

  std::vector<float> max_buffer;
  std::vector<float> buffer;
public:
  Limiter (int n_channels, int sample_rate);

  void set_release (float value_ms);
  void set_attack (float value_ms);
  void set_ceiling (float ceiling);

  std::vector<float> process (const std::vector<float>& samples);
};

#endif /* AUDIOWMARK_LIMITER_HH */
