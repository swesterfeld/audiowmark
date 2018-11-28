#include <string.h>
#include <math.h>
#include <string>
#include <random>
#include <complex>

#include "fft.hh"
#include "wavdata.hh"

using std::string;
using std::vector;
using std::complex;
using std::min;

namespace Params
{
  static constexpr int frame_size      = 1024;
  static constexpr int frames_per_bit  = 4;
  static constexpr int bands_per_frame = 30;
  static constexpr int max_band        = 100;
  static constexpr int min_band        = 20;
  static constexpr double water_gain   = 0.1;   // relative amplitude of the watermark
  static constexpr double pre_scale    = 0.75;  // rescale the signal to avoid clipping after watermark is added
}

inline double
window_cos (double x) /* von Hann window */
{
  if (fabs (x) > 1)
    return 0;
  return 0.5 * cos (x * M_PI) + 0.5;
}

inline double
window_hamming (double x) /* sharp (rectangle) cutoffs at boundaries */
{
  if (fabs (x) > 1)
    return 0;

  return 0.54 + 0.46 * cos (M_PI * x);
}

int
frame_count (WavData& wav_data)
{
  return (wav_data.n_values() / wav_data.n_channels() + (Params::frame_size - 1)) / Params::frame_size;
}

/*
 * get one audio frame, Params::frame_size samples if available
 *
 * in case of stereo: deinterleave
 */
vector<float>
get_frame (WavData& wav_data, int f, int ch)
{
  auto& samples = wav_data.samples();

  vector<float> result;

  size_t pos = (f * Params::frame_size) * wav_data.n_channels() + ch;
  for (size_t x = 0; x < Params::frame_size; x++)
    {
      if (pos < samples.size())
        result.push_back (samples[pos]);

      pos += wav_data.n_channels();
    }
  return result;
}

void
get_up_down (int f, vector<int>& up, vector<int>& down)
{
  vector<int> used (Params::frame_size / 2);
  std::mt19937_64 rng;

  // use per frame random seed, may want to have cryptographically secure algorithm
  rng.seed (f);

  auto choose_bands = [&used, &rng] (vector<int>& bands) {
    while (bands.size() < Params::bands_per_frame)
      {
        int p = rng() % (Params::max_band - Params::min_band) + Params::min_band;
        if (!used[p])
          {
            bands.push_back (p);
            used[p] = 1;
          }
      }
  };
  choose_bands (up);
  choose_bands (down);
}

static unsigned char
from_hex_nibble (char c)
{
  int uc = (unsigned char)c;

  if (uc >= '0' && uc <= '9') return uc - (unsigned char)'0';
  if (uc >= 'a' && uc <= 'f') return uc + 10 - (unsigned char)'a';
  if (uc >= 'A' && uc <= 'F') return uc + 10 - (unsigned char)'A';

  return 16;	// error
}

vector<int>
bit_str_to_vec (const string& bits)
{
  vector<int> bitvec;
  for (auto nibble : bits)
    {
      unsigned char c = from_hex_nibble (nibble);
      if (c >= 16)
        return vector<int>(); // error

      bitvec.push_back ((c & 8) > 0);
      bitvec.push_back ((c & 4) > 0);
      bitvec.push_back ((c & 2) > 0);
      bitvec.push_back ((c & 1) > 0);
    }
  return bitvec;
}

string
bit_vec_to_str (const vector<int>& bit_vec)
{
  string bit_str;

  for (size_t pos = 0; pos < bit_vec.size(); pos += 4)
    {
      int nibble = 0;
      for (int j = 0; j < 4; j++)
        {
          if (pos + j < bit_vec.size())
            if (bit_vec[pos + j])
              {
                // j == 0 has the highest value, then 1, 2, 3 (lowest)
                nibble |= 1 << (3 - j);
              }
        }
      const char *to_hex = "0123456789abcdef";
      bit_str += to_hex[nibble];
    }
  return bit_str;
}

int
add_watermark (const string& infile, const string& outfile, const string& bits)
{
  auto bitvec = bit_str_to_vec (bits);
  if (bitvec.empty())
    {
      fprintf (stderr, "audiowmark: cannot parse bits %s\n", bits.c_str());
      return 1;
    }

  printf ("loading %s\n", infile.c_str());

  WavData wav_data;
  if (!wav_data.load (infile))
    {
      fprintf (stderr, "audiowmark: error loading %s: %s\n", infile.c_str(), wav_data.error_blurb());
      return 1;
    }
  vector<float> out_signal (wav_data.n_values());
  printf ("channels: %d, samples: %zd, mix_freq: %f\n", wav_data.n_channels(), wav_data.n_values(), wav_data.mix_freq());

  vector<float> synth_window (Params::frame_size);
  for (int i = 0; i < Params::frame_size; i++)
    {
      const double threshold = 0.2;

      // triangular basic window
      const double tri = min (1.0 - fabs (double (2 * i)/Params::frame_size - 1.0), threshold) / threshold;

      // cosine
      synth_window[i] = (cos (tri*M_PI+M_PI)+1) * 0.5;

      printf ("cosw %d %f\n", i, synth_window[i]);
    }
  for (int f = 0; f < frame_count (wav_data); f++)
    {
      for (int ch = 0; ch < wav_data.n_channels(); ch++)
        {
          vector<float> frame = get_frame (wav_data, f, ch);
          vector<float> new_frame = frame; /* will be modified below */
          if (frame.size() == Params::frame_size)
            {
              /* windowing */
              double window_weight = 0;
              for (size_t i = 0; i < frame.size(); i++)
                {
                  const double fsize_2 = frame.size() / 2.0;
                  // const double win =  window_cos ((i - fsize_2) / fsize_2);
                  const double win = window_hamming ((i - fsize_2) / fsize_2);
                  //const double win = 1;
                  frame[i] *= win;
                  window_weight += win;
                }

              /* to get normalized fft output corrected by window weight */
              for (size_t i = 0; i < frame.size(); i++)
                frame[i] *= 2.0 / window_weight;

              /* FFT transform */
              vector<complex<float>> fft_out = fft (frame);

              vector<complex<float>> fft_delta_spect (fft_out.size());

              vector<int> up;
              vector<int> down;
              get_up_down (f, up, down);

              const int     data_bit = bitvec[(f / Params::frames_per_bit) % bitvec.size()];
              const double  data_bit_sign = data_bit > 0 ? 1 : -1;
              for (auto u : up)
                {
                  const float factor = Params::water_gain * data_bit_sign;

                  fft_delta_spect[u] = fft_out[u] * factor;
                }
              for (auto d : down)
                {
                  const float factor = -Params::water_gain * data_bit_sign;

                  fft_delta_spect[d] = fft_out[d] * factor;
                }

              for (size_t i = 0; i <= Params::frame_size / 2; i++)
                {
                  const double re = fft_out[i].real();
                  const double im = fft_out[i].imag();
                  const double mag = sqrt (re * re + im * im);
                  printf ("fft %d %d %zd %f\n", f, ch, i, mag);
                }

              vector<float> fft_delta_out = ifft (fft_delta_spect);

              for (size_t i = 0; i < new_frame.size(); i++)
                {
                  printf ("out %d %d %zd %f %f\n", f, ch, i, new_frame[i], fft_delta_out[i]);
                  new_frame[i] += fft_delta_out[i] * synth_window[i];
                }
            }
          for (size_t i = 0; i < new_frame.size(); i++)
            out_signal[(f * Params::frame_size + i) * wav_data.n_channels() + ch] = new_frame[i] * Params::pre_scale;
        }
    }

  WavData out_wav_data (out_signal, wav_data.n_channels(), wav_data.mix_freq(), wav_data.bit_depth());
  if (!out_wav_data.save (outfile))
    {
      fprintf (stderr, "audiowmark: error saving %s: %s\n", outfile.c_str(), wav_data.error_blurb());
      return 1;
    }
  return 0;
}

int
get_watermark (const string& origfile, const string& infile, const string& orig_pattern)
{
  WavData orig_wav_data;
  if (!orig_wav_data.load (origfile))
    {
      fprintf (stderr, "audiowmark: error loading %s: %s\n", origfile.c_str(), orig_wav_data.error_blurb());
      return 1;
    }

  WavData wav_data;
  if (!wav_data.load (infile))
    {
      fprintf (stderr, "audiowmark: error loading %s: %s\n", infile.c_str(), wav_data.error_blurb());
      return 1;
    }
  vector<int> bit_vec;

  double umag = 0, dmag = 0;
  for (int f = 0; f < frame_count (wav_data); f++)
    {
      for (int ch = 0; ch < wav_data.n_channels(); ch++)
        {
          vector<float> frame = get_frame (wav_data, f, ch);
          if (frame.size() == Params::frame_size)
            {
              /* windowing */
              double window_weight = 0;
              for (size_t i = 0; i < frame.size(); i++)
                {
                  const double fsize_2 = frame.size() / 2.0;
                  // const double win =  window_cos ((i - fsize_2) / fsize_2);
                  const double win = window_hamming ((i - fsize_2) / fsize_2);
                  //const double win = 1;
                  frame[i] *= win;
                  window_weight += win;
                }

              /* to get normalized fft output corrected by window weight */
              for (size_t i = 0; i < frame.size(); i++)
                frame[i] *= 2.0 / window_weight;

              /* FFT transform */
              vector<complex<float>> fft_out = fft (frame);

              vector<int> up;
              vector<int> down;
              get_up_down (f, up, down);
              for (auto u : up)
                {
                  umag += log (abs (fft_out[u]));
                }
              for (auto d : down)
                {
                  dmag += log (abs (fft_out[d]));
                }
            }
        }
      if ((f % Params::frames_per_bit) == (Params::frames_per_bit - 1))
        {
          bit_vec.push_back ((umag > dmag) ? 1 : 0);
          umag = 0;
          dmag = 0;
        }
    }
  printf ("pattern %s\n", bit_vec_to_str (bit_vec).c_str());
  if (!orig_pattern.empty())
    {
      int bits = 0, bit_errors = 0;

      vector<int> orig_vec = bit_str_to_vec (orig_pattern);
      for (size_t i = 0; i < bit_vec.size(); i++)
        {
          bits++;
          if (bit_vec[i] != orig_vec[i % orig_vec.size()])
            bit_errors++;
        }
      printf ("bit_error_rate %.3f %%\n", double (100.0 * bit_errors) / bits);
    }
  return 0;
}

int
gentest (const string& infile, const string& outfile)
{
  printf ("generating test sample from '%s' to '%s'\n", infile.c_str(), outfile.c_str());

  WavData wav_data;
  if (!wav_data.load (infile))
    {
      fprintf (stderr, "audiowmark: error loading %s: %s\n", infile.c_str(), wav_data.error_blurb());
      return 1;
    }
  const vector<float>& in_signal = wav_data.samples();
  vector<float> out_signal;

  /* 10 seconds of audio - starting at 30 seconds of the original track */
  const size_t offset = 30 * wav_data.n_channels() * int (wav_data.mix_freq());
  const size_t n_samples = 10 * wav_data.n_channels() * int (wav_data.mix_freq());
  if (in_signal.size() < (offset + n_samples))
    {
      fprintf (stderr, "audiowmark: input file %s too short\n", infile.c_str());
      return 1;
    }
  for (size_t i = 0; i < n_samples; i++)
    {
      out_signal.push_back (in_signal[i + offset]);
    }
  WavData out_wav_data (out_signal, wav_data.n_channels(), wav_data.mix_freq(), wav_data.bit_depth());
  if (!out_wav_data.save (outfile))
    {
      fprintf (stderr, "audiowmark: error saving %s: %s\n", outfile.c_str(), wav_data.error_blurb());
      return 1;
    }
  return 0;
}

int
get_snr (const string& origfile, const string& wmfile)
{
  WavData orig_wav_data;
  if (!orig_wav_data.load (origfile))
    {
      fprintf (stderr, "audiowmark: error loading %s: %s\n", origfile.c_str(), orig_wav_data.error_blurb());
      return 1;
    }

  WavData wav_data;
  if (!wav_data.load (wmfile))
    {
      fprintf (stderr, "audiowmark: error loading %s: %s\n", wmfile.c_str(), wav_data.error_blurb());
      return 1;
    }
  const vector<float>& orig_samples = orig_wav_data.samples();
  const vector<float>& samples = wav_data.samples();

  if (samples.size() != orig_samples.size())
    {
      fprintf (stderr, "audiowmark: files have different length\n");
      return 1;
    }
  double delta_power = 0;
  double signal_power = 0;
  for (size_t i = 0; i < samples.size(); i++)
    {
      const double orig_scaled = orig_samples[i] * Params::pre_scale;
      const double delta       = samples[i] - orig_scaled;

      delta_power += delta * delta;
      signal_power += orig_scaled * orig_scaled;
    }
  delta_power /= samples.size();
  signal_power /= samples.size();

  printf ("snr_db %f\n", 10 * log10 (signal_power / delta_power));
  return 0;
}

int
main (int argc, char **argv)
{
  string op = (argc >= 2) ? argv[1] : "";

  if (op == "add" && argc == 5)
    {
      return add_watermark (argv[2], argv[3], argv[4]);
    }
  else if (op == "get" && argc == 4)
    {
      return get_watermark (argv[2], argv[3], /* no ber */ "");
    }
  else if (op == "cmp" && argc == 5)
    {
      return get_watermark (argv[2], argv[3], argv[4]);
    }
  else if (op == "gentest" && argc == 4)
    {
      return gentest (argv[2], argv[3]);
    }
  else if (op == "snr" && argc == 4)
    {
      get_snr (argv[2], argv[3]);
    }
  else
    {
      fprintf (stderr, "audiowmark: error parsing commandline args\n");
      return 1;
    }
}
