#ifndef AUDIOWMARK_AUDIO_STREAM_HH
#define AUDIOWMARK_AUDIO_STREAM_HH

#include <vector>
#include <memory>
#include "utils.hh"

class AudioStream
{
public:
  virtual int     bit_depth()   const = 0;
  virtual int     sample_rate() const = 0;
  virtual int     n_channels()  const = 0;

  virtual ~AudioStream();
};

class AudioInputStream : public AudioStream
{
public:
  static std::unique_ptr<AudioInputStream> create (const std::string& filename, Error& err);

  // for streams that do not know the number of frames in advance (i.e. raw input stream)
  static constexpr size_t N_FRAMES_UNKNOWN = ~size_t (0);
  virtual size_t n_frames() const = 0;

  virtual Error read_frames (std::vector<float>& samples, size_t count) = 0;
};

class AudioOutputStream : public AudioStream
{
public:
  static std::unique_ptr<AudioOutputStream> create (const std::string& filename,
    int n_channels, int sample_rate, int bit_depth, size_t n_frames, Error& err);

  virtual Error write_frames (const std::vector<float>& frames) = 0;
  virtual Error close() = 0;
};

#endif /* AUDIOWMARK_AUDIO_STREAM_HH */
