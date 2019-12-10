#ifndef AUDIOWMARK_RAW_INPUT_STREAM_HH
#define AUDIOWMARK_RAW_INPUT_STREAM_HH

#include <string>
#include <memory>

#include <sndfile.h>

#include "audiostream.hh"

class RawFormat
{
public:
  enum Endian {
    LITTLE,
    BIG
  };
  enum Encoding {
    SIGNED,
    UNSIGNED
  };
private:
  int       m_n_channels  = 2;
  int       m_sample_rate = 0;
  int       m_bit_depth   = 16;
  Endian    m_endian      = LITTLE;
  Encoding  m_encoding    = SIGNED;
public:
  RawFormat();
  RawFormat (int n_channels, int sample_rate, int bit_depth);

  int n_channels() const { return m_n_channels; }
  int sample_rate() const { return m_sample_rate; }
  int bit_depth() const { return m_bit_depth; }
  Endian endian() const { return m_endian; }
  Encoding encoding() const { return m_encoding; }

  void set_channels (int channels);
  void set_sample_rate (int rate);
  void set_bit_depth (int bits);
  void set_endian (Endian endian);
  void set_encoding (Encoding encoding);
};

class RawConverter;

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

  std::unique_ptr<RawConverter> m_raw_converter;

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
