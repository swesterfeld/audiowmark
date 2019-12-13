#ifndef AUDIOWMARK_STDOUT_WAV_STREAM_HH
#define AUDIOWMARK_STDOUT_WAV_STREAM_HH

#include "audiostream.hh"
#include "rawconverter.hh"

#include <string>

class StdoutWavOutputStream : public AudioOutputStream
{
  int         m_bit_depth = 0;
  int         m_sample_rate = 0;
  int         m_n_channels = 0;
  size_t      m_close_padding = 0;

  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  State       m_state = State::NEW;

  std::unique_ptr<RawConverter> m_raw_converter;

public:
  ~StdoutWavOutputStream();

  Error open (int n_channels, int sample_rate, int bit_depth, size_t n_frames);
  Error write_frames (const std::vector<float>& frames) override;
  Error close();
  int  sample_rate() const override;
  int  bit_depth() const override;
  int  n_channels() const override;
};

#endif
