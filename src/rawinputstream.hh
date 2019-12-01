#ifndef AUDIOWMARK_RAW_INPUT_STREAM_HH
#define AUDIOWMARK_RAW_INPUT_STREAM_HH

#include <string>

#include <sndfile.h>

#include "audiostream.hh"

class RawFormat
{
  int m_n_channels  = 0;
  int m_sample_rate = 0;
  int m_bit_depth   = 0;
public:
  RawFormat();
  RawFormat (int n_channels, int sample_rate, int bit_depth);

  int n_channels() const { return m_n_channels; }
  int sample_rate() const { return m_sample_rate; }
  int bit_depth() const { return m_bit_depth; }

  void set_channels (int channels);
  void set_sample_rate (int rate);
  void set_bit_depth (int bits);
};

class RawInputStream : public AudioInputStream
{
  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  State       m_state = State::NEW;
  RawFormat   m_format;
  FILE       *m_input_file = nullptr;
  bool        m_close_file = false;

public:
  ~RawInputStream();

  Error   open (const std::string& filename, const RawFormat& format);
  Error   read_frames (std::vector<float>& samples, size_t count) override;
  void    close();

  int     bit_depth() const override;
  int     sample_rate() const override;
  size_t  n_frames() const override;
  int     n_channels() const override;
};

#endif /* AUDIOWMARK_RAW_INPUT_STREAM_HH */

