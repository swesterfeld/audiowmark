#ifndef AUDIOWMARK_AUDIO_STREAM_HH
#define AUDIOWMARK_AUDIO_STREAM_HH

#include <vector>

class AudioStream
{
public:
  virtual int     sample_rate() const = 0;
  virtual size_t  n_frames()    const = 0;

  virtual ~AudioStream();
};

class AudioInputStream : public AudioStream
{
  virtual std::vector<float> read_frames (size_t count) = 0;
};

class AudioOutputStream : public AudioStream
{
};

#endif /* AUDIOWMARK_AUDIO_STREAM_HH */
