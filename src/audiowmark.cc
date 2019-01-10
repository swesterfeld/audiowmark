#include <string.h>
#include <math.h>
#include <string>
#include <random>
#include <complex>

#include "fft.hh"
#include "wavdata.hh"
#include "utils.hh"
#include "convcode.hh"
#include "random.hh"

#include <assert.h>

#include "config.h"

using std::string;
using std::vector;
using std::complex;
using std::min;

namespace Params
{
  static size_t frame_size      = 1024;
  static int    frames_per_bit  = 2;
  static size_t bands_per_frame = 30;
  static int max_band          = 100;
  static int min_band          = 20;
  static double water_delta    = 0.015; // strength of the watermark
  static double pre_scale      = 0.95;  // rescale the signal to avoid clipping after watermark is added
  static bool mix              = true;
  static bool hard             = false; // hard decode bits? (soft decoding is better)
  static int block_size        = 32;    // block size for mix step (non-linear bit storage)
  static int have_key          = 0;
  static size_t payload_size   = 128;  // number of payload bits for the watermark

  static int sync_bits           = 6;
  static int sync_frames_per_bit = 32;
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
  printf ("  * generate 128-bit watermarking key, to be used with --key option\n");
  printf ("    audiowmark gen-key <key_file>\n");
  printf ("\n");
  printf ("Global options:\n");
  printf ("  --frame-size          frame size (must be power of 2)     [%zd]\n", Params::frame_size);
  printf ("  --frames-per-bit      number of frames per bit            [%d]\n",  Params::frames_per_bit);
  printf ("  --water-delta         set watermarking delta              [%.4f]\n", Params::water_delta);
  printf ("  --pre-scale           set scaling used for normalization  [%.3f]\n", Params::pre_scale);
  printf ("  --linear              disable non-linear bit storage\n");
  printf ("  --key <file>          load watermarking key from file\n");
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
      else if (check_arg (argc, argv, &i, "--hard"))
	{
          Params::hard = true;
	}
      else if (check_arg (argc, argv, &i, "--test-key", &opt_arg))
	{
          Params::have_key++;
          Random::set_global_test_key (atoi (opt_arg));
	}
      else if (check_arg (argc, argv, &i, "--key", &opt_arg))
        {
          Params::have_key++;
          Random::load_global_key (opt_arg);
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

void
get_up_down (int f, vector<int>& up, vector<int>& down)
{
  vector<int> bands_reorder;
  for (int i = Params::min_band; i <= Params::max_band; i++)
    bands_reorder.push_back (i);

  Random random (f, Random::Stream::up_down); // use per frame random seed
  random.shuffle (bands_reorder);

  assert (2 * Params::bands_per_frame < bands_reorder.size());
  for (size_t i = 0; i < Params::bands_per_frame; i++)
    {
      up.push_back (bands_reorder[i]);
      down.push_back (bands_reorder[Params::bands_per_frame + i]);
    }
}

template<class T> vector<T>
randomize_bit_order (const vector<T>& bit_vec, bool encode)
{
  vector<unsigned int> order;

  for (size_t i = 0; i < bit_vec.size(); i++)
    order.push_back (i);

  Random random (/* seed */ 0, Random::Stream::bit_order);
  random.shuffle (order);

  vector<T> out_bits (bit_vec.size());
  for (size_t i = 0; i < bit_vec.size(); i++)
    {
      if (encode)
        out_bits[i] = bit_vec[order[i]];
      else
        out_bits[order[i]] = bit_vec[i];
    }
  return out_bits;
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
  Random random (/* seed */ block, Random::Stream::mix);
  random.shuffle (mix_entries);

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

size_t
mark_data_frame_count()
{
  const size_t n_blocks = (conv_code_size (Params::payload_size) + (Params::block_size - 1)) / Params::block_size;

  return n_blocks * Params::block_size * Params::frames_per_bit;
}

void
mark_data (const WavData& wav_data, const vector<vector<complex<float>>>& fft_out, vector<vector<complex<float>>>& fft_delta_spect,
           const vector<int>& bitvec)
{
  assert (fft_out.size() == mark_data_frame_count() * wav_data.n_channels());
  assert (bitvec.size() == mark_data_frame_count() / Params::frames_per_bit);

  const int frame_count = fft_out.size() / wav_data.n_channels();

  if (Params::mix)
    {
      const int block_count = frame_count / (Params::block_size * Params::frames_per_bit);

      for (int block = 0; block < block_count; block++)
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

                      const int data_bit = bitvec[(block_start + f) / Params::frames_per_bit];
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
      for (int f = 0; f < frame_count; f++)
        {
          for (int ch = 0; ch < wav_data.n_channels(); ch++)
            {
              size_t index = f * wav_data.n_channels() + ch;

              vector<int> up;
              vector<int> down;
              get_up_down (f, up, down);

              const int data_bit = bitvec[f / Params::frames_per_bit];
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
}

size_t
mark_sync_frame_count()
{
  return Params::sync_bits * Params::sync_frames_per_bit;
}

void
mark_sync (const WavData& wav_data, const vector<vector<complex<float>>>& fft_out, vector<vector<complex<float>>>& fft_delta_spect)
{
  const int frame_count = fft_out.size() / wav_data.n_channels();

  // sync block always written in linear order (no mix)
  for (int f = 0; f < frame_count; f++)
    {
      for (int ch = 0; ch < wav_data.n_channels(); ch++)
        {
          size_t index = f * wav_data.n_channels() + ch;

          vector<int> up;
          vector<int> down;
          get_up_down (f, up, down); // FIXME use extra random stream for sync

          const int data_bit = (f / Params::sync_frames_per_bit) & 1; /* write 010101 */
          const double  data_bit_sign = data_bit > 0 ? 1 : -1;
          for (auto u : up)
            {
              const float mag_factor = pow (abs (fft_out[index][u]), -Params::water_delta * data_bit_sign);

              fft_delta_spect[index][u] = fft_out[index][u] * (mag_factor - 1);
            }
          for (auto d : down)
            {
              const float mag_factor = pow (abs (fft_out[index][d]), Params::water_delta * data_bit_sign);

              fft_delta_spect[index][d] = fft_out[index][d] * (mag_factor - 1);
            }
        }
    }
}

vector<vector<complex<float>>>
get_frame_range (const WavData& wav_data, const vector<vector<complex<float>>>& src, size_t start, size_t count)
{
  start *= wav_data.n_channels();
  count *= wav_data.n_channels();

  assert (start + count < src.size());

  return vector<vector<complex<float>>> (src.begin() + start, src.begin() + start + count);
}

void
copy_frame_range (const WavData& wav_data, const vector<vector<complex<float>>>& src, vector<vector<complex<float>>>& dest, size_t start)
{
  start *= wav_data.n_channels();

  std::copy (src.begin(), src.end(), dest.begin() + start);
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
  if (bitvec.size() > Params::payload_size)
    {
      fprintf (stderr, "audiowmark: number of bits in message '%s' larger than payload size\n", bits.c_str());
      return 1;
    }
  if (bitvec.size() < Params::payload_size)
    {
      /* expand message automatically; good for testing, maybe not so good for the final product */
      vector<int> expanded_bitvec;
      for (size_t i = 0; i < Params::payload_size; i++)
        expanded_bitvec.push_back (bitvec[i % bitvec.size()]);
      bitvec = expanded_bitvec;
    }
  /* add forward error correction, bitvec will now be a lot larger */
  bitvec = randomize_bit_order (conv_encode (bitvec), /* encode */ true);

  /* pad with zeros to match block_size */
  bitvec.resize (mark_data_frame_count() / Params::frames_per_bit);

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
  /*
  vector<vector<complex<float>>> fft_out_range = get_frame_range (wav_data, fft_out, 0, mark_data_frame_count());
  vector<vector<complex<float>>> fft_delta_spect_range = get_frame_range (wav_data, fft_delta_spect, 0, mark_data_frame_count());
  mark_data (wav_data, fft_out_range, fft_delta_spect_range, bitvec);
  copy_frame_range (wav_data, fft_delta_spect_range, fft_delta_spect, 0);
  */
  size_t sync_index = 0;
  while (sync_index + mark_sync_frame_count() < fft_out.size() / wav_data.n_channels())
    {
      vector<vector<complex<float>>> fft_out_range = get_frame_range (wav_data, fft_out, sync_index, mark_sync_frame_count());
      vector<vector<complex<float>>> fft_delta_spect_range = get_frame_range (wav_data, fft_delta_spect, sync_index, mark_sync_frame_count());
      mark_sync (wav_data, fft_out_range, fft_delta_spect_range);
      copy_frame_range (wav_data, fft_delta_spect_range, fft_delta_spect, sync_index);
      sync_index += mark_sync_frame_count();
    }

  /* generate synthesis window */
  // we want overlapping synthesis windows, so the window affects the last, the current and the next frame
  vector<float> synth_window (Params::frame_size * 3);
  for (size_t i = 0; i < synth_window.size(); i++)
    {
      const double overlap = 0.1;

      // triangular basic window
      double tri;
      double norm_pos = (double (i) - Params::frame_size) / Params::frame_size;

      if (norm_pos > 0.5) /* symmetric window */
        norm_pos = 1 - norm_pos;
      if (norm_pos < -overlap)
        {
          tri = 0;
        }
      else if (norm_pos < overlap)
        {
          tri = 0.5 + norm_pos / (2 * overlap);
        }
      else
        {
          tri = 1;
        }
      // cosine
      synth_window[i] = (cos (tri*M_PI+M_PI)+1) * 0.5;
    }

  for (size_t pos = 0; pos < in_signal.size(); pos++)
    out_signal[pos] = in_signal[pos] * Params::pre_scale;

  for (int f = 0; f < frame_count (wav_data); f++)
    {
      for (int ch = 0; ch < wav_data.n_channels(); ch++)
        {
          /* mix watermark signal to output frame */
          vector<float> fft_delta_out = ifft (fft_delta_spect[f * wav_data.n_channels() + ch]);

          int last_frame_start = (f - 1) * Params::frame_size;
          for (int i = 0; i < int (synth_window.size()); i++)
            {
              int pos = (last_frame_start + i) * wav_data.n_channels() + ch;

              if (pos >= 0 && pos < int (out_signal.size()))
                out_signal[pos] += fft_delta_out[i % Params::frame_size] * synth_window[i] * Params::pre_scale;
            }
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

vector<float>
normalize_soft_bits (const vector<float>& soft_bits)
{
  vector<float> norm_soft_bits;

  /* soft decoding produces better error correction than hard decoding */
  if (Params::hard)
    {
      for (auto value : soft_bits)
        norm_soft_bits.push_back (value > 0 ? 1.0 : 0.0);
    }
  else
    {
      /* figure out average level of each bit */
      double mean = 0;
      for (auto value : soft_bits)
        mean += fabs (value);
      mean /= soft_bits.size();

      /* rescale from [-mean,+mean] to [0.0,1.0] */
      for (auto value : soft_bits)
        norm_soft_bits.push_back (0.5 * (value / mean + 1));
    }

  return norm_soft_bits;
}

vector<float>
mix_decode (const WavData& wav_data, vector<vector<complex<float>>>& fft_out, vector<vector<complex<float>>>& fft_orig_out)
{
  vector<float> soft_bit_vec;

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
              soft_bit_vec.push_back (umag - dmag);
              umag = 0;
              dmag = 0;
            }
        }
    }
  return normalize_soft_bits (soft_bit_vec);
}

vector<float>
linear_decode (const WavData& wav_data, vector<vector<complex<float>>>& fft_out, vector<vector<complex<float>>>& fft_orig_out)
{
  vector<float> soft_bit_vec;

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
          soft_bit_vec.push_back (umag - dmag);
          umag = 0;
          dmag = 0;
        }
    }
  return normalize_soft_bits (soft_bit_vec);
}

double
normalize_sync_quality (double raw_quality)
{
  /* the quality for a good sync block depends on watermark strength
   *
   * this is just an approximation, but it should be good enough to be able to
   * use one single threshold on the normalized value check if we have a sync
   * block or not - typical output is 1.0 or more for sync blocks and close
   * to 0.0 for non-sync blocks
   */
  return raw_quality / min (Params::water_delta, 0.080) / 2.9;
}

class SyncFinder
{
  vector<vector<int>> up;
  vector<vector<int>> down;
public:
  SyncFinder()
  {
    up.resize (mark_sync_frame_count());
    down.resize (mark_sync_frame_count());

    for (size_t i = 0; i < mark_sync_frame_count(); i++)
      get_up_down (i, up[i], down[i]);
  }
  double
  sync_decode (const WavData& wav_data, const size_t start_frame, const vector<vector<complex<float>>>& fft_out, const vector<vector<complex<float>>>& fft_orig_out)
  {
    // FIXME: is copypasted
    const int frame_count = mark_sync_frame_count();

    double umag = 0, dmag = 0;
    double sync_quality = 0;

    for (int f = 0; f < frame_count; f++)
      {
        for (int ch = 0; ch < wav_data.n_channels(); ch++)
          {
            const size_t index = (f + start_frame) * wav_data.n_channels() + ch;

            const double min_db = -96;
            for (auto u : up[f])
              {
                umag += db_from_factor (abs (fft_out[index][u]), min_db);

                if (index < fft_orig_out.size())
                  umag -= db_from_factor (abs (fft_orig_out[index][u]), min_db);
              }
            for (auto d : down[f])
              {
                dmag += db_from_factor (abs (fft_out[index][d]), min_db);

                if (index < fft_orig_out.size())
                  dmag -= db_from_factor (abs (fft_orig_out[index][d]), min_db);
              }
          }
        if ((f % Params::sync_frames_per_bit) == (Params::sync_frames_per_bit - 1))
          {
            const int data_bit = (umag < dmag) ? 0 : 1;
            const int expect_data_bit = (f / Params::sync_frames_per_bit) & 1; /* expect 010101 */
            if (data_bit != expect_data_bit)
              return 0;

            const double q = expect_data_bit ? (1 - umag / dmag) : (umag / dmag - 1);
            sync_quality += q;
            umag = 0;
            dmag = 0;
          }
      }
    sync_quality /= Params::sync_bits;
    sync_quality = normalize_sync_quality (sync_quality);
    return sync_quality;
  }
  void
  search (const WavData& wav_data)
  {
    struct SyncScore {
      size_t index;
      double quality;
    };
    vector<SyncScore> sync_scores;

    // compute multiple time-shifted fft vectors
    vector<vector<vector<complex<float>>>> fft_sync_shift_out;
    for (size_t sync_shift = 0; sync_shift < Params::frame_size; sync_shift += 128)
      fft_sync_shift_out.push_back (sync_fft (wav_data, sync_shift, frame_count (wav_data) - 1));

    for (int start_frame = 0; start_frame < frame_count (wav_data); start_frame++)
      {
        for (size_t sync_shift = 0; sync_shift < Params::frame_size; sync_shift += 128)
          {
            const size_t sync_index = start_frame * Params::frame_size + sync_shift;
            if ((start_frame + mark_sync_frame_count()) * wav_data.n_channels() < fft_sync_shift_out[sync_shift / 128].size())
              {
                double quality = sync_decode (wav_data, start_frame, fft_sync_shift_out[sync_shift / 128], /* FIXME: non-blind */ {});
                // printf ("%zd %f\n", sync_index, quality);
                sync_scores.emplace_back (SyncScore { sync_index, quality });
              }
          }
      }
    for (size_t i = 0; i < sync_scores.size(); i++)
      {
        //printf ("%zd %f\n", sync_scores[i].index, sync_scores[i].quality);
        if (sync_scores[i].quality > 0.5)
          {
            double q_last = -1;
            double q_next = -1;

            if (i > 0)
              q_last = sync_scores[i - 1].quality;

            if (i + 1 < sync_scores.size())
              q_next = sync_scores[i + 1].quality;

            if (sync_scores[i].quality > q_last && sync_scores[i].quality > q_next)
              {
                printf ("%zd %s %f", sync_scores[i].index, find_closest_sync (sync_scores[i].index), sync_scores[i].quality);

                // refine match
                double best_quality = sync_scores[i].quality;
                size_t best_index   = sync_scores[i].index;

                int start = std::max (int (sync_scores[i].index) - 128, 0);
                int end   = sync_scores[i].index + 128;
                int step  = 8;
                for (int fine_index = start; fine_index < end; fine_index += step)
                  {
                    vector<vector<complex<float>>> fft_out_range = sync_fft (wav_data, fine_index, mark_sync_frame_count());
                    if (fft_out_range.size())
                      {
                        double q = sync_decode (wav_data, 0, fft_out_range, {});

                        if (q > best_quality)
                          {
                            best_quality = q;
                            best_index   = fine_index;
                          }
                      }
                  }
                printf (" => refined: %zd %s %f\n", best_index, find_closest_sync (best_index), best_quality);
              }
          }
      }
  }
  vector<vector<complex<float>>>
  sync_fft (const WavData& wav_data, size_t index, size_t count)
  {
    if (wav_data.n_values() < (index + count * Params::frame_size) * wav_data.n_channels())
      return {};

    vector<float> part_signal;
    for (size_t i = 0; i < count * Params::frame_size; i++)
      {
        for (int ch = 0; ch < wav_data.n_channels(); ch++)
          part_signal.push_back (wav_data.samples()[(index + i) * wav_data.n_channels() + ch]);
      }
    WavData wav_part (part_signal, wav_data.n_channels(), wav_data.mix_freq(), wav_data.bit_depth());
    return compute_frame_ffts (wav_part);
  }

  const char*
  find_closest_sync (size_t index)
  {
    int best_error = 0xffff;
    int best = 0;

    for (int i = 0; i < 100; i++)
      {
        int error = abs (int (index) - int (i * Params::sync_bits * Params::sync_frames_per_bit * Params::frame_size));
        if (error < best_error)
          {
            best = i;
            best_error = error;
          }
      }
    static char buffer[1024]; // this code is for debugging only, so this should be ok
    sprintf (buffer, "n:%d offset:%d", best, int (index) - int (best * Params::sync_bits * Params::sync_frames_per_bit * Params::frame_size));
    return buffer;
  }
};

int
decode_and_report (const WavData& wav_data, const string& orig_pattern, vector<vector<complex<float>>>& fft_out, vector<vector<complex<float>>>& fft_orig_out)
{
  SyncFinder sync_finder;

  sync_finder.search (wav_data);
  return 0;

  vector<float> soft_bit_vec;
  if (Params::mix)
    {
      soft_bit_vec = mix_decode (wav_data, fft_out, fft_orig_out);
    }
  else
    {
      soft_bit_vec = linear_decode (wav_data, fft_out, fft_orig_out);
    }
  if (soft_bit_vec.size() < conv_code_size (Params::payload_size))
    {
      fprintf (stderr, "audiowmark: input file too short to retrieve watermark\n");
      fprintf (stderr, " - number of recovered raw bits %zd\n", soft_bit_vec.size());
      fprintf (stderr, " - need at least %zd raw bits to get watermark\n", conv_code_size (Params::payload_size));
      return 1;
    }
  /* truncate to the required length */
  soft_bit_vec.resize (conv_code_size (Params::payload_size));

  vector<int> bit_vec = conv_decode_soft (randomize_bit_order (soft_bit_vec, /* encode */ false));

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
      printf ("bit_error_raw %d %d\n", bit_errors, bits);
      printf ("bit_error_rate %.5f %%\n", double (100.0 * bit_errors) / bits);
    }
  return 0;
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

  // to keep the watermark detection code simpler, we truncate samples to avoid partial filled blocks
  truncate_to_block_size (wav_data);

  vector<vector<complex<float>>> fft_out = compute_frame_ffts (wav_data);
  vector<vector<complex<float>>> fft_orig_out; /* no original data -> blind decode */

  return decode_and_report (wav_data, orig_pattern, fft_out, fft_orig_out);
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

  return decode_and_report (wav_data, orig_pattern, fft_out, fft_orig_out);
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

  /* 42 seconds of audio - starting at 30 seconds of the original track */
  /* this is approximately the minimal amount of audio data required for storing a 128-bit encoded message */
  const size_t offset = 30 * wav_data.n_channels() * int (wav_data.mix_freq());
  const size_t n_samples = 42 * wav_data.n_channels() * int (wav_data.mix_freq());
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
cut_start (const string& infile, const string& outfile, const string& start_str)
{
  WavData wav_data;
  if (!wav_data.load (infile))
    {
      fprintf (stderr, "audiowmark: error loading %s: %s\n", infile.c_str(), wav_data.error_blurb());
      return 1;
    }

  size_t start = atoi (start_str.c_str());

  const vector<float>& in_signal = wav_data.samples();
  vector<float> out_signal;
  for (size_t i = start * wav_data.n_channels(); i < in_signal.size(); i++)
    out_signal.push_back (in_signal[i]);

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
gen_key (const string& outfile)
{
  FILE *f = fopen (outfile.c_str(), "w");
  if (!f)
    {
      fprintf (stderr, "audiowmark: error writing to file %s\n", outfile.c_str());
      return 1;
    }
  fprintf (f, "# watermarking key for audiowmark\n\nkey %s\n", Random::gen_key().c_str());
  fclose (f);
  return 0;
}

int
main (int argc, char **argv)
{
  parse_options (&argc, &argv);

  if (Params::have_key > 1)
    {
      fprintf (stderr, "audiowmark: watermark key can at most be set once (--key / --test-key option)\n");
      return 1;
    }
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
  else if (op == "cut-start" && argc == 5)
    {
      cut_start (argv[2], argv[3], argv[4]);
    }
  else if (op == "get-delta" && argc == 4)
    {
      return get_watermark_delta (argv[2], argv[3], /* no ber */ "");
    }
  else if (op == "cmp-delta" && argc == 5)
    {
      return get_watermark_delta (argv[2], argv[3], argv[4]);
    }
  else if (op == "gen-key" && argc == 3)
    {
      return gen_key (argv[2]);
    }
  else
    {
      fprintf (stderr, "audiowmark: error parsing commandline args (use audiowmark -h)\n");
      return 1;
    }
}
