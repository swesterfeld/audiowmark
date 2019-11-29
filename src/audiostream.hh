#ifndef AUDIOWMARK_AUDIO_STREAM_HH
#define AUDIOWMARK_AUDIO_STREAM_HH

#include <vector>
#include "utils.hh"

class AudioStream
{
public:
  // for streams that do not know the number of frames in advance (i.e. raw input stream)
  static constexpr size_t N_FRAMES_UNKNOWN = ~size_t (0);

  virtual int     bit_depth()   const = 0;
  virtual int     sample_rate() const = 0;
  virtual size_t  n_frames()    const = 0;
  virtual int     n_channels()  const = 0;

  virtual ~AudioStream();
};

class AudioInputStream : public AudioStream
{
public:
  virtual Error read_frames (std::vector<float>& samples, size_t count) = 0;
};

class AudioOutputStream : public AudioStream
{
public:
  virtual Error write_frames (const std::vector<float>& frames) = 0;
};

#endif /* AUDIOWMARK_AUDIO_STREAM_HH */
