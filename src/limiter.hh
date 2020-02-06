#ifndef AUDIOWMARK_LIMITER_HH
#define AUDIOWMARK_LIMITER_HH

#include <vector>
#include <sys/types.h>

class Limiter
{
  float maximum       = 1;
  double decay_coeff  = 1;
  uint look_ahead     = 0;

  std::vector<float> max_buffer;
  std::vector<float> buffer;
public:
  Limiter (int sample_rate);

  std::vector<float> process (const std::vector<float>& samples);
};

#endif /* AUDIOWMARK_LIMITER_HH */
