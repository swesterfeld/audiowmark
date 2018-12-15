#include <string.h>
#include <math.h>
#include <string>
#include <random>
#include <complex>

#include "fft.hh"
#include "wavdata.hh"

#include <assert.h>

#include "config.h"

using std::string;
using std::vector;
using std::complex;
using std::min;

namespace Params
{
  static size_t frame_size      = 1024;
  static int    frames_per_bit  = 4;
  static size_t bands_per_frame = 30;
  static int max_band          = 100;
  static int min_band          = 20;
  static double water_delta    = 0.015; // strength of the watermark
  static double pre_scale      = 0.95;  // rescale the signal to avoid clipping after watermark is added
  static bool mix              = true;
  static int block_size        = 32;    // block size for mix step (non-linear bit storage)
  static unsigned int seed     = 0;
}

void
print_usage()
{
  printf ("usage: audiowmark <command> [ <args>... ]\n");
  printf ("\n");
  printf ("Commands:\n");
  printf ("  * create a watermarked wav file with a message\n");
  printf ("    audiowmark add <input_wav> <watermarked_wav> <message_hex>\n");
  printf ("\n");
  printf ("  * blind decoding (retrieve message without original file)\n");
  printf ("    audiowmark get <watermarked_wav>\n");
  printf ("\n");
  printf ("  * compute bit error rate for blind decoding\n");
  printf ("    audiowmark cmp <watermarked_wav> <message_hex>\n");
  printf ("\n");
  printf ("  * retrieve message with original file\n");
  printf ("    audiowmark get-delta <input_wav> <watermarked_wav>\n");
  printf ("\n");
  printf ("  * compute bit error rate for decoding with original file\n");
  printf ("    audiowmark cmp-delta <input_wav> <watermarked_wav> <message_hex>\n");
  printf ("\n");
  printf ("Global options:\n");
  printf ("  --frame-size          frame size (must be power of 2)     [%zd]\n", Params::frame_size);
  printf ("  --frames-per-bit      number of frames per bit            [%d]\n",  Params::frames_per_bit);
  printf ("  --water-delta         set watermarking delta              [%.4f]\n", Params::water_delta);
  printf ("  --pre-scale           set scaling used for normalization  [%.3f]\n", Params::pre_scale);
  printf ("  --linear              disable non-linear bit storage\n");
}

static bool
check_arg (uint         argc,
           char        *argv[],
           uint        *nth,
           const char  *opt,		    /* for example: --foo */
           const char **opt_arg = nullptr)  /* if foo needs an argument, pass a pointer to get the argument */
{
  assert (opt != nullptr);
  assert (*nth < argc);

  const char *arg = argv[*nth];
  if (!arg)
    return false;

  uint opt_len = strlen (opt);
  if (strcmp (arg, opt) == 0)
    {
      if (opt_arg && *nth + 1 < argc)	  /* match foo option with argument: --foo bar */
        {
          argv[(*nth)++] = nullptr;
          *opt_arg = argv[*nth];
          argv[*nth] = nullptr;
          return true;
        }
      else if (!opt_arg)		  /* match foo option without argument: --foo */
        {
          argv[*nth] = nullptr;
          return true;
        }
      /* fall through to error message */
    }
  else if (strncmp (arg, opt, opt_len) == 0 && arg[opt_len] == '=')
    {
      if (opt_arg)			  /* match foo option with argument: --foo=bar */
        {
          *opt_arg = arg + opt_len + 1;
          argv[*nth] = nullptr;
          return true;
        }
      /* fall through to error message */
    }
  else
    return false;

  print_usage();
  exit (1);
}

void
parse_options (int   *argc_p,
               char **argv_p[])
{
  uint argc = *argc_p;
  char **argv = *argv_p;
  unsigned int i, e;

  for (i = 1; i < argc; i++)
    {
      const char *opt_arg;
      if (strcmp (argv[i], "--help") == 0 ||
          strcmp (argv[i], "-h") == 0)
	{
	  print_usage();
	  exit (0);
	}
      else if (strcmp (argv[i], "--version") == 0 || strcmp (argv[i], "-v") == 0)
	{
	  printf ("audiowmark %s\n", VERSION);
	  exit (0);
	}
      else if (check_arg (argc, argv, &i, "--frame-size", &opt_arg))
	{
          Params::frame_size = atoi (opt_arg);
	}
      else if (check_arg (argc, argv, &i, "--frames-per-bit", &opt_arg))
	{
          Params::frames_per_bit = atoi (opt_arg);
	}
      else if (check_arg (argc, argv, &i, "--water-delta", &opt_arg))
	{
          Params::water_delta = atof (opt_arg);
	}
      else if (check_arg (argc, argv, &i, "--pre-scale", &opt_arg))
	{
          Params::pre_scale = atof (opt_arg);
	}
      else if (check_arg (argc, argv, &i, "--linear"))
	{
          Params::mix = false;
	}
      else if (check_arg (argc, argv, &i, "--seed", &opt_arg))
	{
          Params::seed = atoi (opt_arg);
	}
    }

  /* resort argc/argv */
  e = 1;
  for (i = 1; i < argc; i++)
    if (argv[i])
      {
        argv[e++] = argv[i];
        if (i >= e)
          argv[i] = nullptr;
      }
  *argc_p = e;
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

double
db_from_factor (double factor, double min_dB)
{
  if (factor > 0)
    {
      double dB = log10 (factor); /* Bell */
      dB *= 20;
      return dB;
    }
  else
    return min_dB;
}

int
frame_count (const WavData& wav_data)
{
  return (wav_data.n_values() / wav_data.n_channels() + (Params::frame_size - 1)) / Params::frame_size;
}

int
block_count (const WavData& wav_data)
{
  return frame_count (wav_data) / (Params::block_size * Params::frames_per_bit);
}

/*
 * get one audio frame, Params::frame_size samples if available
 *
 * in case of stereo: deinterleave
 */
vector<float>
get_frame (const WavData& wav_data, int f, int ch)
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

std::mt19937_64
init_rng (uint64_t seed)
{
  const uint64_t  prime = 3126986573;

  std::mt19937_64 rng;
  rng.seed (seed + prime * Params::seed);

  return rng;
}

void
get_up_down (int f, vector<int>& up, vector<int>& down)
{
  vector<int> used (Params::frame_size / 2);

  std::mt19937_64 rng = init_rng (f); // use per frame random seed, may want to have cryptographically secure algorithm

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

template<class T> void
gen_shuffle (vector<T>& result, int seed)
{
  std::mt19937_64 rng = init_rng (seed); // should use cryptographically secure generator, properly seeded

  // Fisher–Yates shuffle
  for (size_t i = 0; i < result.size() - 1; i++)
    {
      size_t j = i + rng() % (result.size() - i);
      std::swap (result[i], result[j]);
    }
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

  for (size_t pos = 0; pos + 3 < bit_vec.size(); pos += 4) // convert only groups of 4 bits
    {
      int nibble = 0;
      for (int j = 0; j < 4; j++)
        {
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

struct MixEntry
{
  int  frame;
  int  up;
  int  down;
};

vector<MixEntry>
gen_mix_entries (int block)
{
  vector<MixEntry> mix_entries;

  for (int f = 0; f < Params::block_size * Params::frames_per_bit; f++)
    {
      vector<int> up;
      vector<int> down;
      get_up_down (f, up, down);

      assert (up.size() == down.size());
      for (size_t i = 0; i < up.size(); i++)
        mix_entries.push_back ({ f, up[i], down[i] });
    }
  gen_shuffle (mix_entries, /* seed */ block);
  return mix_entries;
}

vector<vector<complex<float>>>
compute_frame_ffts (const WavData& wav_data)
{
  vector<vector<complex<float>>> fft_out;

  /* generate analysis window */
  vector<float> window (Params::frame_size);

  double window_weight = 0;
  for (size_t i = 0; i < Params::frame_size; i++)
    {
      const double fsize_2 = Params::frame_size / 2.0;
      // const double win =  window_cos ((i - fsize_2) / fsize_2);
      const double win = window_hamming ((i - fsize_2) / fsize_2);
      //const double win = 1;
      window[i] = win;
      window_weight += win;
    }

  /* normalize window using window weight */
  for (size_t i = 0; i < Params::frame_size; i++)
    {
      window[i] *= 2.0 / window_weight;
    }


  for (int f = 0; f < frame_count (wav_data); f++)
    {
      for (int ch = 0; ch < wav_data.n_channels(); ch++)
        {
          vector<float> frame = get_frame (wav_data, f, ch);

          /* apply window */
          for (size_t i = 0; i < frame.size(); i++)
            frame[i] *= window[i];

          /* FFT transform */
          fft_out.push_back (fft (frame));
        }
    }
  return fft_out;
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

  WavData in_wav_data;
  if (!in_wav_data.load (infile))
    {
      fprintf (stderr, "audiowmark: error loading %s: %s\n", infile.c_str(), in_wav_data.error_blurb());
      return 1;
    }

  /*
   * to keep the watermarking code simpler, we pad the wave data with zeros
   * to avoid processing a partly filled block
   */
  vector<float> in_signal (in_wav_data.samples());
  while (in_signal.size() % (in_wav_data.n_channels() * Params::frame_size * Params::block_size * Params::frames_per_bit))
    in_signal.push_back (0);

  WavData wav_data (in_signal, in_wav_data.n_channels(), in_wav_data.mix_freq(), in_wav_data.bit_depth());

  /* we have extra space for the padded wave data -> truncated before save */
  vector<float> out_signal (wav_data.n_values());
  printf ("channels: %d, samples: %zd, mix_freq: %f\n", wav_data.n_channels(), wav_data.n_values(), wav_data.mix_freq());

  vector<vector<complex<float>>> fft_out = compute_frame_ffts (wav_data);
  vector<vector<complex<float>>> fft_delta_spect;
  for (int f = 0; f < frame_count (wav_data); f++)
    {
      for (int ch = 0; ch < wav_data.n_channels(); ch++)
        {
          fft_delta_spect.push_back (vector<complex<float>> (fft_out.back().size()));
        }
    }
  if (Params::mix)
    {
      for (int block = 0; block < block_count (wav_data); block++)
        {
          vector<MixEntry> mix_entries = gen_mix_entries (block);

          const int block_start = block * Params::block_size * Params::frames_per_bit;
          for (int f = 0; f < Params::block_size * Params::frames_per_bit; f++)
            {
              for (int ch = 0; ch < wav_data.n_channels(); ch++)
                {
                  for (size_t frame_b = 0; frame_b < Params::bands_per_frame; frame_b++)
                    {
                      int b = f * Params::bands_per_frame + frame_b;

                      const int data_bit = bitvec[((block_start + f) / Params::frames_per_bit) % bitvec.size()];
                      const double  data_bit_sign = data_bit > 0 ? 1 : -1;

                      const int u = mix_entries[b].up;
                      const int index = (block_start + mix_entries[b].frame) * wav_data.n_channels() + ch;
                      {
                        const float mag_factor = pow (abs (fft_out[index][u]), -Params::water_delta * data_bit_sign);

                        fft_delta_spect[index][u] = fft_out[index][u] * (mag_factor - 1);
                      }
                      const int d = mix_entries[b].down;
                      {
                        const float mag_factor = pow (abs (fft_out[index][d]), Params::water_delta * data_bit_sign);

                        fft_delta_spect[index][d] = fft_out[index][d] * (mag_factor - 1);
                      }
                    }
                }
            }
        }

    }
  else
    {
      for (int f = 0; f < frame_count (wav_data); f++)
        {
          for (int ch = 0; ch < wav_data.n_channels(); ch++)
            {
              size_t index = f * wav_data.n_channels() + ch;

              vector<int> up;
              vector<int> down;
              get_up_down (f, up, down);

              const int data_bit = bitvec[(f / Params::frames_per_bit) % bitvec.size()];
              const double  data_bit_sign = data_bit > 0 ? 1 : -1;
              for (auto u : up)
                {
                  /*
                   * for up bands, we want do use [for a 1 bit]  (pow (mag, 1 - water_delta))
                   *
                   * this actually increases the amount of energy because mag is less than 1.0
                   */
                  const float mag_factor = pow (abs (fft_out[index][u]), -Params::water_delta * data_bit_sign);

                  fft_delta_spect[index][u] = fft_out[index][u] * (mag_factor - 1);
                }
              for (auto d : down)
                {
                  /*
                   * for down bands, we want do use [for a 1 bit]   (pow (mag, 1 + water_delta))
                   *
                   * this actually decreases the amount of energy because mag is less than 1.0
                   */
                  const float mag_factor = pow (abs (fft_out[index][d]), Params::water_delta * data_bit_sign);

                  fft_delta_spect[index][d] = fft_out[index][d] * (mag_factor - 1);
                }
            }
        }
    }

  /* generate synthesis window */
  vector<float> synth_window (Params::frame_size);
  for (size_t i = 0; i < Params::frame_size; i++)
    {
      const double threshold = 0.2;

      // triangular basic window
      const double tri = min (1.0 - fabs (double (2 * i)/Params::frame_size - 1.0), threshold) / threshold;

      // cosine
      synth_window[i] = (cos (tri*M_PI+M_PI)+1) * 0.5;
    }

  for (int f = 0; f < frame_count (wav_data); f++)
    {
      for (int ch = 0; ch < wav_data.n_channels(); ch++)
        {
          /* add watermark to output frame */
          vector<float> frame = get_frame (wav_data, f, ch);

          /* mix watermark signal to output frame */
          vector<float> fft_delta_out = ifft (fft_delta_spect[f * wav_data.n_channels() + ch]);

          for (size_t i = 0; i < frame.size(); i++)
            frame[i] += fft_delta_out[i] * synth_window[i];

          /* modify out signal */
          for (size_t i = 0; i < frame.size(); i++)
            out_signal[(f * Params::frame_size + i) * wav_data.n_channels() + ch] = frame[i] * Params::pre_scale;
        }
    }

  bool clipping_warning = false;
  for (auto value : out_signal)
    {
      if (fabs (value) >= 1.0 && !clipping_warning)
        {
          fprintf (stderr, "audiowmark: warning: clipping occured in watermarked audio signal\n");
          clipping_warning = true;
        }
    }

  out_signal.resize (in_wav_data.n_values()); /* undo zero padding after load */

  WavData out_wav_data (out_signal, wav_data.n_channels(), wav_data.mix_freq(), wav_data.bit_depth());
  if (!out_wav_data.save (outfile))
    {
      fprintf (stderr, "audiowmark: error saving %s: %s\n", outfile.c_str(), out_wav_data.error_blurb());
      return 1;
    }
  return 0;
}

void
truncate_to_block_size (WavData& wav_data)
{
  vector<float> in_signal (wav_data.samples());
  while (in_signal.size() % (wav_data.n_channels() * Params::frame_size * Params::block_size * Params::frames_per_bit))
    in_signal.pop_back();

  wav_data.set_samples (in_signal);
}

vector<int>
mix_decode (const WavData& wav_data, vector<vector<complex<float>>>& fft_out, vector<vector<complex<float>>>& fft_orig_out)
{
  vector<int> bit_vec;

  for (int block = 0; block < block_count (wav_data); block++)
    {
      vector<MixEntry> mix_entries = gen_mix_entries (block);

      double umag = 0, dmag = 0;
      for (int f = 0; f < Params::block_size * Params::frames_per_bit; f++)
        {
          for (int ch = 0; ch < wav_data.n_channels(); ch++)
            {
              for (size_t frame_b = 0; frame_b < Params::bands_per_frame; frame_b++)
                {
                  int b = f * Params::bands_per_frame + frame_b;
                  const double min_db = -96;

                  const size_t index = (block * (Params::block_size * Params::frames_per_bit) + mix_entries[b].frame) * wav_data.n_channels() + ch;
                  const int u = mix_entries[b].up;
                  const int d = mix_entries[b].down;

                  umag += db_from_factor (abs (fft_out[index][u]), min_db);
                  dmag += db_from_factor (abs (fft_out[index][d]), min_db);

                  if (index < fft_orig_out.size()) /* non-blind decode? */
                    {
                      umag -= db_from_factor (abs (fft_orig_out[index][u]), min_db);
                      dmag -= db_from_factor (abs (fft_orig_out[index][d]), min_db);
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
    }
  return bit_vec;
}

vector<int>
linear_decode (const WavData& wav_data, vector<vector<complex<float>>>& fft_out, vector<vector<complex<float>>>& fft_orig_out)
{
  vector<int> bit_vec;

  double umag = 0, dmag = 0;
  for (int f = 0; f < frame_count (wav_data); f++)
    {
      for (int ch = 0; ch < wav_data.n_channels(); ch++)
        {
          const size_t index = f * wav_data.n_channels() + ch;
          vector<int> up;
          vector<int> down;
          get_up_down (f, up, down);

          const double min_db = -96;
          for (auto u : up)
            {
              umag += db_from_factor (abs (fft_out[index][u]), min_db);

              if (index < fft_orig_out.size())
                umag -= db_from_factor (abs (fft_orig_out[index][u]), min_db);
            }
          for (auto d : down)
            {
              dmag += db_from_factor (abs (fft_out[index][d]), min_db);

              if (index < fft_orig_out.size())
                dmag -= db_from_factor (abs (fft_orig_out[index][d]), min_db);
            }
        }
      if ((f % Params::frames_per_bit) == (Params::frames_per_bit - 1))
        {
          bit_vec.push_back ((umag > dmag) ? 1 : 0);
          umag = 0;
          dmag = 0;
        }
    }
  return bit_vec;
}

int
get_watermark (const string& infile, const string& orig_pattern)
{
  WavData wav_data;
  if (!wav_data.load (infile))
    {
      fprintf (stderr, "audiowmark: error loading %s: %s\n", infile.c_str(), wav_data.error_blurb());
      return 1;
    }
  vector<int> bit_vec;

  // to keep the watermark detection code simpler, we truncate samples to avoid partial filled blocks
  truncate_to_block_size (wav_data);

  vector<vector<complex<float>>> fft_out = compute_frame_ffts (wav_data);
  vector<vector<complex<float>>> fft_orig_out; /* no original data -> blind decode */

  if (Params::mix)
    {
      bit_vec = mix_decode (wav_data, fft_out, fft_orig_out);
    }
  else
    {
      bit_vec = linear_decode (wav_data, fft_out, fft_orig_out);
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
      printf ("bit_error_rate %.5f %%\n", double (100.0 * bit_errors) / bits);
    }
  return 0;
}

int
get_watermark_delta (const string& origfile, const string& infile, const string& orig_pattern)
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

  // to keep the watermark detection code simpler, we truncate samples to avoid partial filled blocks
  truncate_to_block_size (wav_data);
  truncate_to_block_size (orig_wav_data);

  vector<vector<complex<float>>> fft_out = compute_frame_ffts (wav_data);
  vector<vector<complex<float>>> fft_orig_out = compute_frame_ffts (orig_wav_data);

  vector<int> bit_vec;
  if (Params::mix)
    {
      bit_vec = mix_decode (wav_data, fft_out, fft_orig_out);
    }
  else
    {
      bit_vec = linear_decode (wav_data, fft_out, fft_orig_out);
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
      printf ("bit_error_rate %.5f %%\n", double (100.0 * bit_errors) / bits);
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
      fprintf (stderr, "audiowmark: error saving %s: %s\n", outfile.c_str(), out_wav_data.error_blurb());
      return 1;
    }
  return 0;
}

int
scale (const string& infile, const string& outfile)
{
  WavData wav_data;
  if (!wav_data.load (infile))
    {
      fprintf (stderr, "audiowmark: error loading %s: %s\n", infile.c_str(), wav_data.error_blurb());
      return 1;
    }

  const vector<float>& in_signal = wav_data.samples();
  vector<float> out_signal;
  for (size_t i = 0; i < in_signal.size(); i++)
    out_signal.push_back (in_signal[i] * Params::pre_scale);

  WavData out_wav_data (out_signal, wav_data.n_channels(), wav_data.mix_freq(), wav_data.bit_depth());
  if (!out_wav_data.save (outfile))
    {
      fprintf (stderr, "audiowmark: error saving %s: %s\n", outfile.c_str(), out_wav_data.error_blurb());
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
  parse_options (&argc, &argv);

  string op = (argc >= 2) ? argv[1] : "";

  if (op == "add" && argc == 5)
    {
      return add_watermark (argv[2], argv[3], argv[4]);
    }
  else if (op == "get" && argc == 3)
    {
      return get_watermark (argv[2], /* no ber */ "");
    }
  else if (op == "cmp" && argc == 4)
    {
      return get_watermark (argv[2], argv[3]);
    }
  else if (op == "gentest" && argc == 4)
    {
      return gentest (argv[2], argv[3]);
    }
  else if (op == "snr" && argc == 4)
    {
      get_snr (argv[2], argv[3]);
    }
  else if (op == "scale" && argc == 4)
    {
      scale (argv[2], argv[3]);
    }
  else if (op == "get-delta" && argc == 4)
    {
      return get_watermark_delta (argv[2], argv[3], /* no ber */ "");
    }
  else if (op == "cmp-delta" && argc == 5)
    {
      return get_watermark_delta (argv[2], argv[3], argv[4]);
    }
  else
    {
      fprintf (stderr, "audiowmark: error parsing commandline args\n");
      return 1;
    }
}
