#ifndef AUDIOWMARK_WAV_DATA_HH
#define AUDIOWMARK_WAV_DATA_HH

#include <string>
#include <vector>

#include "utils.hh"
#include "audiostream.hh"

class WavData
{
  std::vector<float> m_samples;
  int                m_sample_rate = 0;
  int                m_n_channels  = 0;
  int                m_bit_depth   = 0;

public:
  WavData();
  WavData (const std::vector<float>& samples, int n_channels, int sample_rate, int bit_depth);

  Error load (AudioInputStream *in_stream);
  Error load (const std::string& filename);
  Error save (const std::string& filename);

  int                         sample_rate() const;
  int                         bit_depth() const;

  int
  n_channels() const
  {
    return m_n_channels;
  }
  size_t
  n_values() const
  {
    return m_samples.size();
  }
  const std::vector<float>&
  samples() const
  {
    return m_samples;
  }

  void set_samples (const std::vector<float>& samples);
};

#endif /* AUDIOWMARK_WAV_DATA_HH */
