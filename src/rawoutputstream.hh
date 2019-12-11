#ifndef AUDIOWMARK_RAW_OUTPUT_STREAM_HH
#define AUDIOWMARK_RAW_OUTPUT_STREAM_HH

#include "rawinputstream.hh"
#include "rawconverter.hh"

#include <memory>

class RawOutputStream : public AudioOutputStream
{
  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  State       m_state = State::NEW;
  RawFormat   m_format;
  FILE       *m_output_file = nullptr;
  bool        m_close_file = false;

  std::unique_ptr<RawConverter> m_raw_converter;
public:
  ~RawOutputStream();

  int   bit_depth() const override;
  int   sample_rate() const override;
  int   n_channels()  const override;

  Error open (const std::string& filename, const RawFormat& format);
  Error write_frames (const std::vector<float>& frames) override;
  Error close();
};

#endif /* AUDIOWMARK_RAW_OUTPUT_STREAM_HH */
