#ifndef AUDIOWMARK_LIMITER_HH
#define AUDIOWMARK_LIMITER_HH

#include <vector>
#include <sys/types.h>

class Limiter
{
  double ceiling        = 1;
  double maximum        = 1;
  double release_factor = 0;
  uint look_ahead       = 0;
  uint sample_rate      = 0;

  std::vector<float> max_buffer;
  std::vector<float> buffer;
public:
  Limiter (int sample_rate);

  void set_release (double value_ms);
  void set_attack (double value_ms);
  void set_ceiling (double ceiling);

  std::vector<float> process (const std::vector<float>& samples);
};

#endif /* AUDIOWMARK_LIMITER_HH */
