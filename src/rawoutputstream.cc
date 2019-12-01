#include "rawoutputstream.hh"

#include <assert.h>

using std::string;
using std::vector;

RawOutputStream::~RawOutputStream()
{
  close();
}

Error
RawOutputStream::open (const string& filename, const RawFormat& format)
{
  assert (m_state == State::NEW);

  if (!format.n_channels())
    return Error ("RawOutputStream: output format: missing number of channels");
  if (!format.bit_depth())
    return Error ("RawOutputStream: output format: missing bit depth");
  if (!format.sample_rate())
    return Error ("RawOutputStream: output format: missing sample rate");

  Error err = Error::Code::NONE;
  RawConverter *rc = RawConverter::create (format, err);
  if (err)
    return err;
  assert (rc);
  m_raw_converter.reset (rc);

  if (filename == "-")
    {
      m_output_file = stdout;
      m_close_file = false;
    }
  else
    {
      m_output_file = fopen (filename.c_str(), "w");
      m_close_file = true;
    }

  m_format = format;
  m_state  = State::OPEN;
  return Error::Code::NONE;
}

int
RawOutputStream::sample_rate() const
{
  return m_format.sample_rate();
}

int
RawOutputStream::bit_depth() const
{
  return m_format.bit_depth();
}

int
RawOutputStream::n_channels() const
{
  return m_format.n_channels();
}

Error
RawOutputStream::write_frames (const vector<float>& samples)
{
  assert (m_state == State::OPEN);

  vector<unsigned char> bytes;
  m_raw_converter->to_raw (samples, bytes);

  fwrite (&bytes[0], 1, bytes.size(), m_output_file);

  return Error::Code::NONE;
}

void
RawOutputStream::close()
{
  if (m_state == State::OPEN)
    {
      if (m_close_file && m_output_file)
        {
          fclose (m_output_file);
          m_output_file = nullptr;
          m_close_file = false;
        }

      m_state = State::CLOSED;
    }
}
