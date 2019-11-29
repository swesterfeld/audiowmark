#ifndef AUDIOWMARK_RAW_INPUT_STREAM_HH
#define AUDIOWMARK_RAW_INPUT_STREAM_HH

#include <string>

#include <sndfile.h>

#include "audiostream.hh"

class RawInputStream : public AudioInputStream
{
  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  State       m_state = State::NEW;

public:
  ~RawInputStream();

  Error   open (const std::string& filename);
  Error   read_frames (std::vector<float>& samples, size_t count) override;
  void    close();

  int     bit_depth() const override;
  int     sample_rate() const override;
  size_t  n_frames() const override;
  int     n_channels() const override;
};

#endif /* AUDIOWMARK_RAW_INPUT_STREAM_HH */

