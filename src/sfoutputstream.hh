#ifndef AUDIOWMARK_SF_OUTPUT_STREAM_HH
#define AUDIOWMARK_SF_OUTPUT_STREAM_HH

#include <string>

#include <sndfile.h>

#include "audiostream.hh"

class SFOutputStream : public AudioOutputStream
{
  SNDFILE    *m_sndfile = nullptr;
  int         m_bit_depth = 0;
  int         m_sample_rate = 0;
  int         m_n_channels = 0;

  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  State       m_state = State::NEW;

public:
  ~SFOutputStream();

  Error  open (const std::string& filename, int n_channels, int sample_rate, int bit_depth, size_t n_frames);
  Error  write_frames (const std::vector<float>& frames) override;
  void   close();
  int    bit_depth() const override;
  int    sample_rate() const override;
  int    n_channels() const override;
};

#endif /* AUDIOWMARK_SF_OUTPUT_STREAM_HH */
