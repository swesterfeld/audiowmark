#include <string>
#include <vector>
#include <sndfile.h>
#include <assert.h>
#include <math.h>

#include "wavdata.hh"
#include "utils.hh"

class AudioInputStream
{
  virtual std::vector<float> read_frames (size_t count) = 0;
};

class AudioOutputStream
{
};

class SFInputStream : public AudioInputStream
{
  SNDFILE    *m_sndfile = nullptr;
  std::string m_error_blurb;
  int         m_n_channels = 0;
  int         m_n_values = 0;
  int         m_sample_rate = 0;

  enum class State {
    NEW,
    OPEN,
    CLOSED
  };
  State       m_state = State::NEW;

public:
  ~SFInputStream();

  bool                open (const std::string& filename);
  std::vector<float>  read_frames (size_t count);
  void                close();

  int
  n_channels() const
  {
    return m_n_channels;
  }
  int sample_rate() const;
  size_t
  n_values() const
  {
    return m_n_values;
  }
  size_t
  n_frames() const
  {
    return m_n_values / m_n_channels;
  }
  const char *error_blurb() const
  {
    return m_error_blurb.c_str();
  }
};

class StdoutWavOutputStream : public AudioOutputStream
{
  std::string m_error_blurb;
public:
  bool open (int n_channels, int sample_rate, int bit_depth, size_t n_frames);
  bool write_frames (const std::vector<float>& frames);

  const char *error_blurb() const
  {
    return m_error_blurb.c_str();
  }
};

using std::string;
using std::vector;

SFInputStream::~SFInputStream()
{
  close();
}

bool
SFInputStream::open (const string& filename)
{
  assert (m_state == State::NEW);

  SF_INFO sfinfo = { 0, };

  m_sndfile = sf_open (filename.c_str(), SFM_READ, &sfinfo);

  int error = sf_error (m_sndfile);
  if (error)
    {
      m_error_blurb = sf_strerror (m_sndfile);
      if (m_sndfile)
        {
          m_sndfile = nullptr;
          sf_close (m_sndfile);
        }
      return false;
    }

  m_n_channels  = sfinfo.channels;
  m_n_values    = sfinfo.frames * sfinfo.channels;
  m_sample_rate = sfinfo.samplerate;
  m_state       = State::OPEN;
  return true;
}

int
SFInputStream::sample_rate() const
{
  return m_sample_rate;
}

vector<float>
SFInputStream::read_frames (size_t count)
{
  assert (m_state == State::OPEN);

  vector<int> isamples (count * m_n_channels);
  sf_count_t r_count = sf_readf_int (m_sndfile, &isamples[0], count);

  /* reading a wav file and saving it again with the libsndfile float API will
   * change some values due to normalization issues:
   *   http://www.mega-nerd.com/libsndfile/FAQ.html#Q010
   *
   * to avoid the problem, we use the int API and do the conversion beween int
   * and float manually - the important part is that the normalization factors
   * used during read and write are identical
   */
  vector<float> result (r_count * m_n_channels);;
  const double norm = 1.0 / 0x80000000LL;
  for (size_t i = 0; i < result.size(); i++)
    result[i] = isamples[i] * norm;

  return result;
}

void
SFInputStream::close()
{
  if (m_state == State::OPEN)
    {
      assert (m_sndfile);
      sf_close (m_sndfile);

      m_state = State::CLOSED;
    }
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

bool
StdoutWavOutputStream::open (int n_channels, int sample_rate, int bit_depth, size_t n_frames)
{
  if (bit_depth != 16)
    {
      m_error_blurb = "StdoutWavOutputStream::open: unsupported bit depth";
      return false;
    }
  vector<unsigned char> header_bytes;

  size_t data_size = n_frames * n_channels * ((bit_depth + 7) / 8);

  header_append_str (header_bytes, "RIFF");
  header_append_u32 (header_bytes, 36 + data_size);
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
  return true;
}

bool
StdoutWavOutputStream::write_frames (const vector<float>& samples)
{
  vector<unsigned char> output_bytes (samples.size() * sizeof (short));

  for (size_t i = 0; i < samples.size(); i++)
    {
      const double norm      =  0x80000000LL;
      const double min_value = -0x80000000LL;
      const double max_value =  0x7FFFFFFF;

      const int    sample = lrint (bound<double> (min_value, samples[i] * norm, max_value));

      // write short little endian value
      output_bytes[i * 2]     = sample >> 16;
      output_bytes[i * 2 + 1] = sample >> 24;
    }
  fwrite (&output_bytes[0], 1, output_bytes.size(), stdout);
  return true;
}

int
main (int argc, char **argv)
{
  SFInputStream in;
  StdoutWavOutputStream out;

  std::string filename = (argc >= 2) ? argv[1] : "-";
  if (!in.open (filename.c_str()))
    {
      fprintf (stderr, "teststream: open input failed: %s\n", in.error_blurb());
      return 1;
    }
  if (!out.open (in.n_channels(), in.sample_rate(), 16, in.n_frames()))
    {
      fprintf (stderr, "teststream: open output failed: %s\n", out.error_blurb());
      return 1;
    }
  vector<float> samples;
  do
    {
      samples = in.read_frames (1024);
      out.write_frames (samples);
    }
  while (samples.size());
}
