#ifndef AUDIOWMARK_MP3_INPUT_STREAM_HH
#define AUDIOWMARK_MP3_INPUT_STREAM_HH

#include <string>
#include <mpg123.h>

#include "audiostream.hh"

class MP3InputStream : public AudioInputStream
{
  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  int         m_n_values = 0;
  int         m_n_channels = 0;
  int         m_sample_rate = 0;
  int         m_frames_left = 0;
  bool        m_need_close = false;
  State       m_state = State::NEW;

  mpg123_handle     *m_handle = nullptr;
  std::vector<float> m_read_buffer;
public:
  ~MP3InputStream();

  Error   open (const std::string& filename);
  Error   read_frames (std::vector<float>& samples, size_t count);
  void    close();

  int     bit_depth() const override;
  int     sample_rate() const override;
  int     n_channels()  const override;
  size_t  n_frames() const override;

};

#endif /* AUDIOWMARK_MP3_INPUT_STREAM_HH */

