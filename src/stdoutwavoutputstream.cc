#include "stdoutwavoutputstream.hh"
#include "utils.hh"

#include <assert.h>
#include <math.h>

using std::string;
using std::vector;

StdoutWavOutputStream::~StdoutWavOutputStream()
{
  close();
}

int
StdoutWavOutputStream::sample_rate() const
{
  return m_sample_rate;
}

int
StdoutWavOutputStream::bit_depth() const
{
  return m_bit_depth;
}

static void
header_append_str (vector<unsigned char>& bytes, const string& str)
{
  for (auto ch : str)
    bytes.push_back (ch);
}

static void
header_append_u32 (vector<unsigned char>& bytes, uint32_t u)
{
  bytes.push_back (u);
  bytes.push_back (u >> 8);
  bytes.push_back (u >> 16);
  bytes.push_back (u >> 24);
}

static void
header_append_u16 (vector<unsigned char>& bytes, uint16_t u)
{
  bytes.push_back (u);
  bytes.push_back (u >> 8);
}

Error
StdoutWavOutputStream::open (int n_channels, int sample_rate, int bit_depth, size_t n_frames)
{
  assert (m_state == State::NEW);

  if (bit_depth != 16 && bit_depth != 24)
    {
      return Error ("StdoutWavOutputStream::open: unsupported bit depth");
    }
  vector<unsigned char> header_bytes;

  size_t data_size = n_frames * n_channels * ((bit_depth + 7) / 8);

  m_close_padding = data_size & 1; // padding to ensure even data size
  size_t aligned_data_size = data_size + m_close_padding;

  header_append_str (header_bytes, "RIFF");
  header_append_u32 (header_bytes, 36 + aligned_data_size);
  header_append_str (header_bytes, "WAVE");

  // subchunk 1
  header_append_str (header_bytes, "fmt ");
  header_append_u32 (header_bytes, 16); // subchunk size
  header_append_u16 (header_bytes, 1);  // uncompressed audio
  header_append_u16 (header_bytes, n_channels);
  header_append_u32 (header_bytes, sample_rate);
  header_append_u32 (header_bytes, sample_rate * n_channels * bit_depth / 8); // byte rate
  header_append_u16 (header_bytes, n_channels * bit_depth / 8); // block align
  header_append_u16 (header_bytes, bit_depth); // bits per sample

  // subchunk 2
  header_append_str (header_bytes, "data");
  header_append_u32 (header_bytes, data_size);

  fwrite (&header_bytes[0], 1, header_bytes.size(), stdout);

  m_bit_depth = bit_depth;
  m_n_frames  = n_frames;
  m_state     = State::OPEN;

  return Error::Code::NONE;
}

template<int BIT_DEPTH> void
convert_frames (const vector<float>& samples, vector<unsigned char>& output_bytes)
{
  for (size_t i = 0; i < samples.size(); i++)
    {
      const double norm      =  0x80000000LL;
      const double min_value = -0x80000000LL;
      const double max_value =  0x7FFFFFFF;

      const int    sample = lrint (bound<double> (min_value, samples[i] * norm, max_value));

      if (BIT_DEPTH == 16)
        {
          // write 16-bit little endian value
          output_bytes[i * 2]     = sample >> 16;
          output_bytes[i * 2 + 1] = sample >> 24;
        }
      else if (BIT_DEPTH == 24)
        {
          // write 24-bit little endian value
          output_bytes[i * 3]     = sample >> 8;
          output_bytes[i * 3 + 1] = sample >> 16;
          output_bytes[i * 3 + 2] = sample >> 24;
        }
    }
}

Error
StdoutWavOutputStream::write_frames (const vector<float>& samples)
{
  vector<unsigned char> output_bytes (samples.size() * (m_bit_depth / 8));

  if (m_bit_depth == 16)
    {
      convert_frames<16> (samples, output_bytes);
    }
  else if (m_bit_depth == 24)
    {
      convert_frames<24> (samples, output_bytes);
    }
  else
    {
      assert (false);
    }

  fwrite (&output_bytes[0], 1, output_bytes.size(), stdout);
  return Error::Code::NONE;
}

void
StdoutWavOutputStream::close()
{
  if (m_state == State::OPEN)
    {
      for (size_t i = 0; i < m_close_padding; i++)
        fputc (0, stdout);

      m_state = State::CLOSED;
    }
}


