#ifndef AUDIOWMARK_SF_INPUT_STREAM_HH
#define AUDIOWMARK_SF_INPUT_STREAM_HH

#include <string>

#include <sndfile.h>

#include "audiostream.hh"

class SFInputStream : public AudioInputStream
{
  SNDFILE    *m_sndfile = nullptr;
  std::string m_error_blurb;
  int         m_n_channels = 0;
  int         m_n_values = 0;
  int         m_bit_depth = 0;
  int         m_sample_rate = 0;

  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  State       m_state = State::NEW;

public:
  ~SFInputStream();

  bool                open (const std::string& filename);
  Error               read_frames (std::vector<float>& samples, size_t count);
  void                close();

  int
  n_channels() const
  {
    return m_n_channels;
  }
  int sample_rate() const override;
  int bit_depth() const override;
  size_t
  n_values() const
  {
    return m_n_values;
  }
  size_t
  n_frames() const override
  {
    return m_n_values / m_n_channels;
  }
  const char *error_blurb() const
  {
    return m_error_blurb.c_str();
  }
};

#endif /* AUDIOWMARK_SF_INPUT_STREAM_HH */
