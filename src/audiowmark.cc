#include <string.h>
#include <math.h>
#include <string>
#include <random>
#include <complex>
#include <algorithm>

#include "fft.hh"
#include "wavdata.hh"
#include "utils.hh"
#include "convcode.hh"
#include "random.hh"

#include <zita-resampler/resampler.h>
#include <zita-resampler/vresampler.h>

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
  static double water_delta    = 0.01;  // strength of the watermark
  static double pre_scale      = 0.95;  // rescale the signal to avoid clipping after watermark is added
  static bool mix              = true;
  static bool hard             = false; // hard decode bits? (soft decoding is better)
  static int have_key          = 0;
  static size_t payload_size   = 128;  // number of payload bits for the watermark

  static int sync_bits           = 6;
  static int sync_frames_per_bit = 85;
  static int sync_search_step    = 256;
  static int sync_search_fine    = 8;
  static double sync_threshold1  = 0.5; // minimum grid quality value (search_step grid)
  static double sync_threshold2  = 0.7; // minimum refined quality

  static size_t frames_pad_start = 250; // padding at start, in case track starts with silence
  static int    mark_sample_rate = 44100; // watermark generation and detection sample rate

  static int test_cut            = 0; // for sync test
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
#if 0
  printf ("  * retrieve message with original file\n");
  printf ("    audiowmark get-delta <input_wav> <watermarked_wav>\n");
  printf ("\n");
  printf ("  * compute bit error rate for decoding with original file\n");
  printf ("    audiowmark cmp-delta <input_wav> <watermarked_wav> <message_hex>\n");
  printf ("\n");
#endif
  printf ("  * generate 128-bit watermarking key, to be used with --key option\n");
  printf ("    audiowmark gen-key <key_file>\n");
  printf ("\n");
  printf ("Global options:\n");
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
      else if (check_arg (argc, argv, &i, "--test-cut", &opt_arg))
	{
          Params::test_cut = atoi (opt_arg);
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
  return wav_data.n_values() / wav_data.n_channels() / Params::frame_size;
}

void
get_up_down (int f, vector<int>& up, vector<int>& down, Random::Stream random_stream)
{
  vector<int> bands_reorder;
  for (int i = Params::min_band; i <= Params::max_band; i++)
    bands_reorder.push_back (i);

  Random random (f, random_stream); // use per frame random seed
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

vector<vector<complex<float>>>
compute_frame_ffts (const WavData& wav_data, size_t start_index, size_t frame_count)
{
  vector<vector<complex<float>>> fft_out;

  /* if there is not enough space for frame_count values, return an error (empty vector) */
  if (wav_data.n_values() < (start_index + frame_count * Params::frame_size) * wav_data.n_channels())
    return fft_out;

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

  float *frame  = new_array_float (Params::frame_size);
  float *frame_fft = new_array_float (Params::frame_size);

  for (size_t f = 0; f < frame_count; f++)
    {
      for (int ch = 0; ch < wav_data.n_channels(); ch++)
        {
          const auto& samples = wav_data.samples();

          size_t pos = (start_index + f * Params::frame_size) * wav_data.n_channels() + ch;
          assert (pos + (Params::frame_size - 1) * wav_data.n_channels() < samples.size());

          /* deinterleave frame data and apply window */
          for (size_t x = 0; x < Params::frame_size; x++)
            {
              frame[x] = samples[pos] * window[x];
              pos += wav_data.n_channels();
            }
          /* FFT transform */
          fftar_float (Params::frame_size, frame, frame_fft);

          /* complex<float> and frame_fft have the same layout in memory */
          const complex<float> *first = (complex<float> *) frame_fft;
          const complex<float> *last  = first + Params::frame_size / 2 + 1;
          fft_out.emplace_back (first, last);
        }
    }
  free_array_float (frame);
  free_array_float (frame_fft);
  return fft_out;
}

void
mark_bit_linear (int f, const vector<complex<float>>& fft_out, vector<complex<float>>& fft_delta_spect, int data_bit, Random::Stream random_stream)
{
  vector<int> up;
  vector<int> down;
  get_up_down (f, up, down, random_stream);
  const double  data_bit_sign = data_bit > 0 ? 1 : -1;
  for (auto u : up)
    {
      /*
       * for up bands, we want do use [for a 1 bit]  (pow (mag, 1 - water_delta))
       *
       * this actually increases the amount of energy because mag is less than 1.0
       */
      const float mag_factor = pow (abs (fft_out[u]), -Params::water_delta * data_bit_sign);

      fft_delta_spect[u] = fft_out[u] * (mag_factor - 1);
    }
  for (auto d : down)
    {
      /*
       * for down bands, we want do use [for a 1 bit]   (pow (mag, 1 + water_delta))
       *
       * this actually decreases the amount of energy because mag is less than 1.0
       */
      const float mag_factor = pow (abs (fft_out[d]), Params::water_delta * data_bit_sign);

      fft_delta_spect[d] = fft_out[d] * (mag_factor - 1);
    }
}

size_t
mark_data_frame_count()
{
  return conv_code_size (Params::payload_size) * Params::frames_per_bit;
}

struct MixEntry
{
  int  frame;
  int  up;
  int  down;
};

vector<MixEntry>
gen_mix_entries()
{
  vector<MixEntry> mix_entries;

  for (int f = 0; f < int (mark_data_frame_count()); f++)
    {
      vector<int> up;
      vector<int> down;
      get_up_down (f, up, down, Random::Stream::data_up_down);

      assert (up.size() == down.size());
      for (size_t i = 0; i < up.size(); i++)
        mix_entries.push_back ({ f, up[i], down[i] });
    }
  Random random (/* seed */ 0, Random::Stream::mix);
  random.shuffle (mix_entries);

  return mix_entries;
}

void
mark_data (const WavData& wav_data, int start_frame, const vector<vector<complex<float>>>& fft_out, vector<vector<complex<float>>>& fft_delta_spect,
           const vector<int>& bitvec)
{
  assert (fft_out.size() >= (start_frame + mark_data_frame_count()) * wav_data.n_channels());
  assert (bitvec.size() == mark_data_frame_count() / Params::frames_per_bit);

  const int frame_count = mark_data_frame_count();

  if (Params::mix)
    {
      vector<MixEntry> mix_entries = gen_mix_entries();

      for (int f = 0; f < frame_count; f++)
        {
          for (int ch = 0; ch < wav_data.n_channels(); ch++)
            {
              for (size_t frame_b = 0; frame_b < Params::bands_per_frame; frame_b++)
                {
                  int b = f * Params::bands_per_frame + frame_b;

                  const int data_bit = bitvec[f / Params::frames_per_bit];
                  const double  data_bit_sign = data_bit > 0 ? 1 : -1;

                  const int u = mix_entries[b].up;
                  const int index = (start_frame + mix_entries[b].frame) * wav_data.n_channels() + ch;
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
  else
    {
      for (int f = 0; f < frame_count; f++)
        {
          for (int ch = 0; ch < wav_data.n_channels(); ch++)
            {
              size_t index = (start_frame + f) * wav_data.n_channels() + ch;

              mark_bit_linear (f, fft_out[index], fft_delta_spect[index], bitvec[f / Params::frames_per_bit], Random::Stream::data_up_down);
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
mark_sync (const WavData& wav_data, int start_frame, const vector<vector<complex<float>>>& fft_out, vector<vector<complex<float>>>& fft_delta_spect)
{
  assert (fft_out.size() >= (start_frame + mark_sync_frame_count()) * wav_data.n_channels());

  const int frame_count = mark_sync_frame_count();

  // sync block always written in linear order (no mix)
  for (int f = 0; f < frame_count; f++)
    {
      for (int ch = 0; ch < wav_data.n_channels(); ch++)
        {
          size_t index = (start_frame + f) * wav_data.n_channels() + ch;
          int    data_bit = (f / Params::sync_frames_per_bit) & 1; /* write 010101 */

          mark_bit_linear (f, fft_out[index], fft_delta_spect[index], data_bit, Random::Stream::sync_up_down);
        }
    }
}

void
mark_pad (const WavData& wav_data, size_t frame, const vector<vector<complex<float>>>& fft_out, vector<vector<complex<float>>>& fft_delta_spect)
{
  assert (fft_out.size() >= (frame + 1) * wav_data.n_channels());

  for (int ch = 0; ch < wav_data.n_channels(); ch++)
    {
      size_t index = frame * wav_data.n_channels() + ch;

      mark_bit_linear (frame, fft_out[index], fft_delta_spect[index], 0, Random::Stream::pad_up_down);
    }
}

template<class R>
static void
process_resampler (R& resampler, const vector<float>& in, vector<float>& out)
{
  resampler.out_count = out.size() / resampler.nchan();
  resampler.out_data = &out[0];

  /* avoid timeshift: zita needs k/2 - 1 samples before the actual input */
  resampler.inp_count = resampler.inpsize () / 2 - 1;
  resampler.inp_data  = nullptr;
  resampler.process();

  resampler.inp_count = in.size() / resampler.nchan();
  resampler.inp_data = (float *) &in[0];
  resampler.process();

  /* zita needs k/2 samples after the actual input */
  resampler.inp_count = resampler.inpsize() / 2;
  resampler.inp_data  = nullptr;
  resampler.process();
}

WavData
resample (const WavData& wav_data, int rate)
{
  /* in our application, resampling should only be called if it is necessary
   * since using the resampler with input rate == output rate would be slow
   */
  assert (rate != wav_data.sample_rate());

  const int hlen = 16;
  const double ratio = double (rate) / wav_data.sample_rate();

  const vector<float>& in = wav_data.samples();
  vector<float> out (lrint (in.size() / wav_data.n_channels() * ratio) * wav_data.n_channels());

  /* zita-resampler provides two resampling algorithms
   *
   * a fast optimized version: Resampler
   *   this is an optimized version, which works for many common cases,
   *   like resampling between 22050, 32000, 44100, 48000, 96000 Hz
   *
   * a slower version: VResampler
   *   this works for arbitary rates (like 33333 -> 44100 resampling)
   *
   * so we try using Resampler, and if that fails fall back to VResampler
   */
  Resampler resampler;
  if (resampler.setup (wav_data.sample_rate(), rate, wav_data.n_channels(), hlen) == 0)
    {
      process_resampler (resampler, in, out);
      return WavData (out, wav_data.n_channels(), rate, wav_data.bit_depth());
    }

  VResampler vresampler;
  if (vresampler.setup (ratio, wav_data.n_channels(), hlen) == 0)
    {
      process_resampler (vresampler, in, out);
      return WavData (out, wav_data.n_channels(), rate, wav_data.bit_depth());
    }
  fprintf (stderr, "audiowmark: resampling from rate %d to rate %d not supported.\n", wav_data.sample_rate(), rate);
  exit (1);
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
  printf ("Input:        %s\n", infile.c_str());
  printf ("Output:       %s\n", outfile.c_str());
  printf ("Message:      %s\n\n", bit_vec_to_str (bitvec).c_str());

  /* add forward error correction, bitvec will now be a lot larger */
  bitvec = randomize_bit_order (conv_encode (bitvec), /* encode */ true);

  /* pad with zeros to match block_size */
  bitvec.resize (mark_data_frame_count() / Params::frames_per_bit);

  WavData orig_wav_data;
  if (!orig_wav_data.load (infile))
    {
      fprintf (stderr, "audiowmark: error loading %s: %s\n", infile.c_str(), orig_wav_data.error_blurb());
      return 1;
    }
  int orig_seconds = orig_wav_data.n_values() / orig_wav_data.sample_rate() / orig_wav_data.n_channels();
  printf ("Time:         %d:%02d\n", orig_seconds / 60, orig_seconds % 60);
  printf ("Sample Rate:  %d\n", orig_wav_data.sample_rate());
  printf ("Channels:     %d\n", orig_wav_data.n_channels());

  vector<float> in_signal;
  if (orig_wav_data.sample_rate() != Params::mark_sample_rate)
    {
      WavData in_wav_data = resample (orig_wav_data, Params::mark_sample_rate);
      in_signal = in_wav_data.samples();
    }
  else
    {
      in_signal = orig_wav_data.samples();
    }

  /*
   * to keep the watermarking code simpler, we pad the wave data with zeros
   * to avoid processing a partly filled frame
   */
  while (in_signal.size() % (orig_wav_data.n_channels() * Params::frame_size))
    in_signal.push_back (0);

  WavData wav_data (in_signal, orig_wav_data.n_channels(), Params::mark_sample_rate, orig_wav_data.bit_depth());

  /* we have extra space for the padded wave data -> truncated before save */
  vector<float> out_signal (wav_data.n_values());

  vector<vector<complex<float>>> fft_out = compute_frame_ffts (wav_data, 0, frame_count (wav_data));
  vector<vector<complex<float>>> fft_delta_spect;
  for (int f = 0; f < frame_count (wav_data); f++)
    {
      for (int ch = 0; ch < wav_data.n_channels(); ch++)
        {
          fft_delta_spect.push_back (vector<complex<float>> (fft_out.back().size()));
        }
    }
  size_t frame_index = 0;
  int    data_blocks = 0;
  /* padding at start */
  while (frame_index < Params::frames_pad_start)
    {
      mark_pad (wav_data, frame_index, fft_out, fft_delta_spect);
      frame_index++;
    }
  /* embed sync|data|sync|data|... */
  while (frame_index + (mark_sync_frame_count() + mark_data_frame_count()) < size_t (frame_count (wav_data)))
    {
      mark_sync (wav_data, frame_index, fft_out, fft_delta_spect);
      frame_index += mark_sync_frame_count();

      data_blocks++;

      mark_data (wav_data, frame_index, fft_out, fft_delta_spect, bitvec);
      frame_index += mark_data_frame_count();
    }
  /* padding at end */
  while (frame_index < size_t (frame_count (wav_data)))
    {
      mark_pad (wav_data, frame_index, fft_out, fft_delta_spect);
      frame_index++;
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
  for (int f = 0; f < frame_count (wav_data); f++)
    {
      for (int ch = 0; ch < wav_data.n_channels(); ch++)
        {
          /* mix watermark signal to output frame */
          vector<float> fft_delta_out = ifft (fft_delta_spect[f * wav_data.n_channels() + ch]);

          for (int dframe = -1; dframe <= 1; dframe++)
            {
              const int wstart = (dframe + 1) * Params::frame_size;

              if (f + dframe > 0 && f + dframe < frame_count (wav_data))
                {
                  int pos = (f + dframe) * Params::frame_size * wav_data.n_channels() + ch;

                  for (size_t x = 0; x < Params::frame_size; x++)
                    {
                      out_signal[pos] += fft_delta_out[x] * synth_window[wstart + x];
                      pos += wav_data.n_channels();
                    }
                }
            }
        }
    }

  if (wav_data.sample_rate() != orig_wav_data.sample_rate())
    {
      /* resample the watermark to the original sample rate */
      WavData mark_wav_data (out_signal, wav_data.n_channels(), wav_data.sample_rate(), wav_data.bit_depth());
      mark_wav_data = resample (mark_wav_data, orig_wav_data.sample_rate());

      out_signal = mark_wav_data.samples();
    }
  vector<float> samples = orig_wav_data.samples();
  out_signal.resize (samples.size());

  for (size_t i = 0; i < samples.size(); i++)
    samples[i] = (samples[i] + out_signal[i]) * Params::pre_scale;

  bool clipping_warning = false;
  for (auto value : samples)
    {
      if (fabs (value) >= 1.0 && !clipping_warning)
        {
          fprintf (stderr, "audiowmark: warning: clipping occured in watermarked audio signal\n");
          clipping_warning = true;
        }
    }
  printf ("Data Blocks:  %d\n", data_blocks);

  WavData out_wav_data (samples, orig_wav_data.n_channels(), orig_wav_data.sample_rate(), orig_wav_data.bit_depth());
  if (!out_wav_data.save (outfile))
    {
      fprintf (stderr, "audiowmark: error saving %s: %s\n", outfile.c_str(), out_wav_data.error_blurb());
      return 1;
    }
  return 0;
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
mix_decode (vector<vector<complex<float>>>& fft_out, vector<vector<complex<float>>>& fft_orig_out, int n_channels)
{
  vector<float> raw_bit_vec;

  const int frame_count = fft_out.size() / n_channels;

  vector<MixEntry> mix_entries = gen_mix_entries();

  double umag = 0, dmag = 0;
  for (int f = 0; f < frame_count; f++)
    {
      for (int ch = 0; ch < n_channels; ch++)
        {
          for (size_t frame_b = 0; frame_b < Params::bands_per_frame; frame_b++)
            {
              int b = f * Params::bands_per_frame + frame_b;
              const double min_db = -96;

              const size_t index = mix_entries[b].frame * n_channels + ch;
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
          raw_bit_vec.push_back (umag - dmag);
          umag = 0;
          dmag = 0;
        }
    }
  return raw_bit_vec;
}

vector<float>
linear_decode (vector<vector<complex<float>>>& fft_out, vector<vector<complex<float>>>& fft_orig_out, int n_channels)
{
  vector<float> raw_bit_vec;

  double umag = 0, dmag = 0;
  const int frame_count = fft_out.size() / n_channels;
  for (int f = 0; f < frame_count; f++)
    {
      for (int ch = 0; ch < n_channels; ch++)
        {
          const size_t index = f * n_channels + ch;
          vector<int> up;
          vector<int> down;
          get_up_down (f, up, down, Random::Stream::data_up_down);

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
          raw_bit_vec.push_back (umag - dmag);
          umag = 0;
          dmag = 0;
        }
    }
  return raw_bit_vec;
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

  void
  init_up_down (const WavData& wav_data)
  {
    up.clear();
    down.clear();

    up.resize (Params::sync_bits);
    down.resize (Params::sync_bits);

    size_t n_bands = Params::max_band - Params::min_band + 1;
    for (int bit = 0; bit < Params::sync_bits; bit++)
      {
        for (int f = 0; f < Params::sync_frames_per_bit; f++)
          {
            vector<int> frame_up, frame_down;
            get_up_down (f + bit * Params::sync_frames_per_bit, frame_up, frame_down, Random::Stream::sync_up_down);

            for (auto u : frame_up)
              up[bit].push_back (u - Params::min_band + f * n_bands * wav_data.n_channels());

            for (auto d : frame_down)
              down[bit].push_back (d - Params::min_band + f * n_bands * wav_data.n_channels());
          }
        sort (up[bit].begin(), up[bit].end());
        sort (down[bit].begin(), down[bit].end());
      }
  }
  double
  sync_decode (const WavData& wav_data, const size_t start_frame, const vector<float>& fft_out_db)
  {
    double sync_quality = 0;

    size_t n_bands = Params::max_band - Params::min_band + 1;
    for (int bit = 0; bit < Params::sync_bits; bit++)
      {
        float umag = 0, dmag = 0;

        for (int ch = 0; ch < wav_data.n_channels(); ch++)
          {
            const int index = (((bit * Params::sync_frames_per_bit) + start_frame) * wav_data.n_channels() + ch) * n_bands;

            for (size_t i = 0; i < up[bit].size(); i++)
              {
                umag += fft_out_db[index + up[bit][i]];
                dmag += fft_out_db[index + down[bit][i]];
              }
          }
        /* convert avoiding bias, raw_bit < 0 => 0 bit received; raw_bit > 0 => 1 bit received */
        double raw_bit;
        if (umag < dmag)
          {
            raw_bit = 1 - umag / dmag;
          }
        else
          {
            raw_bit = dmag / umag - 1;
          }

        const int expect_data_bit = bit & 1; /* expect 010101 */
        const double q = expect_data_bit ? raw_bit : -raw_bit;
        sync_quality += q;
      }
    sync_quality /= Params::sync_bits;
    sync_quality = normalize_sync_quality (sync_quality);
    return sync_quality;
  }
public:
  struct Score {
    size_t index;
    double quality;
  };
  vector<Score>
  search (const WavData& wav_data)
  {
    vector<Score> result_scores;
    vector<Score> sync_scores;

    init_up_down (wav_data);

    vector<float> fft_db;

    // compute multiple time-shifted fft vectors
    size_t n_bands = Params::max_band - Params::min_band + 1;
    for (size_t sync_shift = 0; sync_shift < Params::frame_size; sync_shift += Params::sync_search_step)
      {
        sync_fft (wav_data, sync_shift, frame_count (wav_data) - 1, fft_db);
        for (int start_frame = 0; start_frame < frame_count (wav_data); start_frame++)
          {
            const size_t sync_index = start_frame * Params::frame_size + sync_shift;
            if ((start_frame + mark_sync_frame_count()) * wav_data.n_channels() * n_bands < fft_db.size())
              {
                double quality = sync_decode (wav_data, start_frame, fft_db);
                // printf ("%zd %f\n", sync_index, quality);
                sync_scores.emplace_back (Score { sync_index, quality });
              }
          }
      }
    sort (sync_scores.begin(), sync_scores.end(), [] (const Score& a, const Score &b) { return a.index < b.index; });
    for (size_t i = 0; i < sync_scores.size(); i++)
      {
        // printf ("%zd %f\n", sync_scores[i].index, sync_scores[i].quality);
        if (sync_scores[i].quality > Params::sync_threshold1)
          {
            double q_last = -1;
            double q_next = -1;

            if (i > 0)
              q_last = sync_scores[i - 1].quality;

            if (i + 1 < sync_scores.size())
              q_next = sync_scores[i + 1].quality;

            if (sync_scores[i].quality > q_last && sync_scores[i].quality > q_next)
              {
                //printf ("%zd %s %f", sync_scores[i].index, find_closest_sync (sync_scores[i].index), sync_scores[i].quality);

                // refine match
                double best_quality = sync_scores[i].quality;
                size_t best_index   = sync_scores[i].index;

                int start = std::max (int (sync_scores[i].index) - Params::sync_search_step, 0);
                int end   = sync_scores[i].index + Params::sync_search_step;
                for (int fine_index = start; fine_index <= end; fine_index += Params::sync_search_fine)
                  {
                    sync_fft (wav_data, fine_index, mark_sync_frame_count(), fft_db);
                    if (fft_db.size())
                      {
                        double q = sync_decode (wav_data, 0, fft_db);

                        if (q > best_quality)
                          {
                            best_quality = q;
                            best_index   = fine_index;
                          }
                      }
                  }
                //printf (" => refined: %zd %s %f\n", best_index, find_closest_sync (best_index), best_quality);
                if (best_quality > Params::sync_threshold2)
                  result_scores.push_back (Score { best_index, best_quality });
              }
          }
      }
    return result_scores;
  }
private:
  void
  sync_fft (const WavData& wav_data, size_t index, size_t count, vector<float>& fft_out_db)
  {
    fft_out_db.clear();

    /* computing db-magnitude is expensive, so we better do it here */
    for (const vector<complex<float>>& spect : compute_frame_ffts (wav_data, index, count))
      {
        const double min_db = -96;

        for (int i = Params::min_band; i <= Params::max_band; i++)
          fft_out_db.push_back (db_from_factor (abs (spect[i]), min_db));
      }
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
decode_and_report (const WavData& wav_data, const string& orig_pattern)
{
  int match_count = 0, total_count = 0, sync_match = 0;

  SyncFinder                sync_finder;
  vector<SyncFinder::Score> sync_scores = sync_finder.search (wav_data);

  auto decode_single = [&] (const vector<float>& raw_bit_vec, SyncFinder::Score sync_score)
  {
    assert (raw_bit_vec.size() == conv_code_size (Params::payload_size));

    vector<float> soft_bit_vec = normalize_soft_bits (raw_bit_vec);

    float decode_error = 0;
    vector<int> bit_vec = conv_decode_soft (randomize_bit_order (soft_bit_vec, /* encode */ false), &decode_error);

    if (sync_score.index)
      {
        const int seconds = sync_score.index / wav_data.sample_rate();
        printf ("pattern %2d:%02d %s %.3f %.3f\n", seconds / 60, seconds % 60, bit_vec_to_str (bit_vec).c_str(), sync_score.quality, decode_error);
      }
    else /* this is the combined pattern "all" */
      {
        printf ("pattern   all %s %.3f %.3f\n", bit_vec_to_str (bit_vec).c_str(), sync_score.quality, decode_error);
      }

    if (!orig_pattern.empty())
      {
        bool        match = true;
        vector<int> orig_vec = bit_str_to_vec (orig_pattern);

        for (size_t i = 0; i < bit_vec.size(); i++)
          match = match && (bit_vec[i] == orig_vec[i % orig_vec.size()]);

        if (match)
          match_count++;

      }
    total_count++;
  };

  vector<float> raw_bit_vec_all (conv_code_size (Params::payload_size));
  SyncFinder::Score score_all { 0, 0 };

  for (auto sync_score : sync_scores)
    {
      const size_t count = mark_data_frame_count();
      const size_t index = sync_score.index + (mark_sync_frame_count() * Params::frame_size);

      auto fft_range_out = compute_frame_ffts (wav_data, index, count);
      if (fft_range_out.size())
        {
          vector<vector<complex<float>>> junk;

          vector<float> raw_bit_vec;
          if (Params::mix)
            {
              raw_bit_vec = mix_decode (fft_range_out, junk, wav_data.n_channels());
            }
          else
            {
              raw_bit_vec = linear_decode (fft_range_out, junk, wav_data.n_channels());
            }
          decode_single (raw_bit_vec, sync_score);

          score_all.quality += sync_score.quality;
          for (size_t i = 0; i < raw_bit_vec_all.size(); i++)
            raw_bit_vec_all[i] += raw_bit_vec[i];
        }
    }
  if (total_count > 1) /* all pattern: average soft bits of all watermarks and decode */
    {
      score_all.quality /= total_count;
      decode_single (raw_bit_vec_all, score_all);
    }

  if (!orig_pattern.empty())
    {
      printf ("match_count %d %d\n", match_count, total_count);

      /* search sync markers at typical positions */
      const int expect0 = Params::frames_pad_start * Params::frame_size;
      const int expect_step = (mark_sync_frame_count() + mark_data_frame_count()) * Params::frame_size;
      const int expect_end = frame_count (wav_data) * Params::frame_size;

      for (int expect_index = expect0; expect_index + expect_step < expect_end; expect_index += expect_step)
        {
          for (auto sync_score : sync_scores)
            {
              if (abs (int (sync_score.index + Params::test_cut) - expect_index) < Params::frame_size / 2)
                {
                  sync_match++;
                  break;
                }
            }
        }
      printf ("sync_match %d %zd\n", sync_match, sync_scores.size());
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

  if (wav_data.sample_rate() == Params::mark_sample_rate)
    {
      return decode_and_report (wav_data, orig_pattern);
    }
  else
    {
      return decode_and_report (resample (wav_data, Params::mark_sample_rate), orig_pattern);
    }
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

  /*
  vector<vector<complex<float>>> fft_out = compute_frame_ffts (wav_data);
  vector<vector<complex<float>>> fft_orig_out = compute_frame_ffts (orig_wav_data);

  return decode_and_report (wav_data, orig_pattern, fft_out, fft_orig_out);
  */
  /* FIXME? */
  printf ("delta decoding currently not supported\n");
  return 1;
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

  /* 2:45 of audio - this is approximately the minimal amount of audio data required
   * for storing three separate watermarks with a 128-bit encoded message */
  const size_t offset = 0 * wav_data.n_channels() * wav_data.sample_rate();
  const size_t n_samples = 165 * wav_data.n_channels() * wav_data.sample_rate();
  if (in_signal.size() < (offset + n_samples))
    {
      fprintf (stderr, "audiowmark: input file %s too short\n", infile.c_str());
      return 1;
    }
  for (size_t i = 0; i < n_samples; i++)
    {
      out_signal.push_back (in_signal[i + offset]);
    }
  WavData out_wav_data (out_signal, wav_data.n_channels(), wav_data.sample_rate(), wav_data.bit_depth());
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

  WavData out_wav_data (out_signal, wav_data.n_channels(), wav_data.sample_rate(), wav_data.bit_depth());
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

  WavData out_wav_data (out_signal, wav_data.n_channels(), wav_data.sample_rate(), wav_data.bit_depth());
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
