#ifndef AUDIOWMARK_STDOUT_WAV_STREAM_HH
#define AUDIOWMARK_STDOUT_WAV_STREAM_HH

#include "audiostream.hh"

#include <string>

class StdoutWavOutputStream : public AudioOutputStream
{
  std::string m_error_blurb;
  int         m_bit_depth = 0;
  int         m_sample_rate = 0;
  size_t      m_n_frames = 0;
  size_t      m_close_padding = 0;

  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  State       m_state = State::NEW;

public:
  ~StdoutWavOutputStream();

  bool open (int n_channels, int sample_rate, int bit_depth, size_t n_frames);
  bool write_frames (const std::vector<float>& frames) override;
  void close();
  int  sample_rate() const override;
  int  bit_depth() const override;
  size_t n_frames() const override
  {
    return m_n_frames;
  }

  const char *error_blurb() const
  {
    return m_error_blurb.c_str();
  }
};

#endif
