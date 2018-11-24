#ifndef AUDIOWMARK_WAV_DATA_HH
#define AUDIOWMARK_WAV_DATA_HH

#include <string>
#include <vector>

class WavData
{
  std::vector<float> m_samples;
  float              m_mix_freq   = 0;
  int                m_n_channels = 0;
  int                m_bit_depth  = 0;
  std::string        m_error_blurb;

public:
  bool load (const std::string& filename);
  bool save (const std::string& filename);

  float                       mix_freq() const;
  int                         n_channels() const;
  size_t                      n_values() const;
  int                         bit_depth() const;
  const std::vector<float>&   samples() const;
  const char                 *error_blurb() const;
};

#endif /* AUDIOWMARK_WAV_DATA_HH */
