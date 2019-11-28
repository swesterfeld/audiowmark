#include <string.h>
#include <math.h>
#include <string>
#include <random>
#include <complex>
#include <algorithm>
#include <memory>

#include "fft.hh"
#include "wavdata.hh"
#include "utils.hh"
#include "convcode.hh"
#include "random.hh"
#include "sfinputstream.hh"
#include "sfoutputstream.hh"
#include "stdoutwavoutputstream.hh"

#include <zita-resampler/resampler.h>
#include <zita-resampler/vresampler.h>

#include <assert.h>

#include "config.h"

using std::string;
using std::vector;
using std::complex;
using std::min;
using std::max;

namespace Params
{
  static size_t frame_size      = 1024;
  static int    frames_per_bit  = 2;
  static size_t bands_per_frame = 30;
  static int max_band          = 100;
  static int min_band          = 20;
  static double water_delta    = 0.01;  // strength of the watermark
  static bool mix              = true;
  static bool hard             = false; // hard decode bits? (soft decoding is better)
  static bool snr              = false; // compute/show snr while adding watermark
  static int have_key          = 0;
  static size_t payload_size   = 128;  // number of payload bits for the watermark

  static int sync_bits           = 6;
  static int sync_frames_per_bit = 85;
  static int sync_search_step    = 256;
  static int sync_search_fine    = 8;
  static double sync_threshold2  = 0.7; // minimum refined quality

  static size_t frames_pad_start = 250; // padding at start, in case track starts with silence
  static int    mark_sample_rate = 44100; // watermark generation and detection sample rate

  static int test_cut            = 0; // for sync test
  static bool test_no_sync       = false; // disable sync
  static int test_truncate       = 0;
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
  printf ("  * retrieve message\n");
  printf ("    audiowmark get <watermarked_wav>\n");
  printf ("\n");
  printf ("  * compare watermark message with expected message\n");
  printf ("    audiowmark cmp <watermarked_wav> <message_hex>\n");
  printf ("\n");
  printf ("  * generate 128-bit watermarking key, to be used with --key option\n");
  printf ("    audiowmark gen-key <key_file>\n");
  printf ("\n");
  printf ("Global options:\n");
  printf ("  --strength <s>        set watermark strength              [%.6g]\n", Params::water_delta * 1000);
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
      else if (check_arg (argc, argv, &i, "--strength", &opt_arg))
	{
          Params::water_delta = atof (opt_arg) / 1000;
	}
      else if (check_arg (argc, argv, &i, "--linear"))
	{
          Params::mix = false;
	}
      else if (check_arg (argc, argv, &i, "--hard"))
	{
          Params::hard = true;
	}
      else if (check_arg (argc, argv, &i, "--snr"))
        {
          Params::snr = true;
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
      else if (check_arg (argc, argv, &i, "--test-no-sync"))
        {
          Params::test_no_sync = true;
        }
      else if (check_arg (argc, argv, &i, "--test-truncate", &opt_arg))
	{
          Params::test_truncate = atoi (opt_arg);
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

/*
 * glibc log2f is a lot faster than glibc log10
 */
inline double
log10_approx (double l)
{
  constexpr double log2_log10_factor = 0.3010299956639811952; // 1 / log2 (10)

  return log2f (l) * log2_log10_factor;
}

double
db_from_factor (double factor, double min_dB)
{
  if (factor > 0)
    {
      double dB = log10_approx (factor); /* Bell */
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

typedef std::array<int, 30> UpDownArray;
class UpDownGen
{
  Random::Stream random_stream;
  Random         random;
  vector<int>    bands_reorder;

public:
  UpDownGen (Random::Stream random_stream) :
    random_stream (random_stream),
    random (0, random_stream),
    bands_reorder (Params::max_band - Params::min_band + 1)
  {
    UpDownArray x;
    assert (x.size() == Params::bands_per_frame);
  }
  void
  get (int f, UpDownArray& up, UpDownArray& down)
  {
    for (size_t i = 0; i < bands_reorder.size(); i++)
      bands_reorder[i] = Params::min_band + i;

    random.seed (f, random_stream); // use per frame random seed
    random.shuffle (bands_reorder);

    assert (2 * Params::bands_per_frame < bands_reorder.size());
    for (size_t i = 0; i < Params::bands_per_frame; i++)
      {
        up[i]   = bands_reorder[i];
        down[i] = bands_reorder[Params::bands_per_frame + i];
      }
  }
};

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

class FFTAnalyzer
{
  int           m_n_channels = 0;
  vector<float> m_window;
  float        *m_frame = nullptr;
  float        *m_frame_fft = nullptr;
public:
  FFTAnalyzer (int n_channels) :
    m_n_channels (n_channels)
  {
    /* generate analysis window */
    m_window.resize (Params::frame_size);

    double window_weight = 0;
    for (size_t i = 0; i < Params::frame_size; i++)
      {
        const double fsize_2 = Params::frame_size / 2.0;
        // const double win =  window_cos ((i - fsize_2) / fsize_2);
        const double win = window_hamming ((i - fsize_2) / fsize_2);
        //const double win = 1;
        m_window[i] = win;
        window_weight += win;
      }

    /* normalize window using window weight */
    for (size_t i = 0; i < Params::frame_size; i++)
      {
        m_window[i] *= 2.0 / window_weight;
      }

    /* allocate properly aligned buffers for SIMD */
    m_frame  = new_array_float (Params::frame_size);
    m_frame_fft = new_array_float (Params::frame_size);
  }
  ~FFTAnalyzer()
  {
    free_array_float (m_frame);
    free_array_float (m_frame_fft);
  }
  vector<vector<complex<float>>>
  run_fft (const vector<float>& samples, size_t start_index)
  {
    assert (samples.size() >= (Params::frame_size + start_index) * m_n_channels);

    vector<vector<complex<float>>> fft_out;
    for (int ch = 0; ch < m_n_channels; ch++)
      {
        size_t pos = start_index * m_n_channels + ch;
        assert (pos + (Params::frame_size - 1) * m_n_channels < samples.size());

        /* deinterleave frame data and apply window */
        for (size_t x = 0; x < Params::frame_size; x++)
          {
            m_frame[x] = samples[pos] * m_window[x];
            pos += m_n_channels;
          }
        /* FFT transform */
        fftar_float (Params::frame_size, m_frame, m_frame_fft);

        /* complex<float> and frame_fft have the same layout in memory */
        const complex<float> *first = (complex<float> *) m_frame_fft;
        const complex<float> *last  = first + Params::frame_size / 2 + 1;
        fft_out.emplace_back (first, last);
      }

    return fft_out;
  }
  vector<vector<complex<float>>>
  fft_range (const vector<float>& samples, size_t start_index, size_t frame_count)
  {
    vector<vector<complex<float>>> fft_out;

    /* if there is not enough space for frame_count values, return an error (empty vector) */
    if (samples.size() < (start_index + frame_count * Params::frame_size) * m_n_channels)
      return fft_out;

    for (size_t f = 0; f < frame_count; f++)
      {
        const size_t frame_start = (f * Params::frame_size) + start_index;

        vector<vector<complex<float>>> frame_result = run_fft (samples, frame_start);
        for (auto& fr : frame_result)
          fft_out.emplace_back (std::move (fr));
      }
    return fft_out;
  }
};

size_t mark_data_frame_count();
size_t mark_sync_frame_count();

int
frame_pos (int f, bool sync)
{
  static vector<int> pos_vec;

  if (pos_vec.empty())
    {
      int frame_count = mark_data_frame_count() + mark_sync_frame_count();
      for (int i = 0; i < frame_count; i++)
        pos_vec.push_back (i);

      Random random (0, Random::Stream::frame_position);
      random.shuffle (pos_vec);
    }
  if (sync)
    {
      assert (f >= 0 && size_t (f) < mark_sync_frame_count());

      return pos_vec[f];
    }
  else
    {
      assert (f >= 0 && size_t (f) < mark_data_frame_count());

      return pos_vec[f + mark_sync_frame_count()];
    }
}

int
sync_frame_pos (int f)
{
  return frame_pos (f, true);
}

int
data_frame_pos (int f)
{
  return frame_pos (f, false);
}

size_t
mark_data_frame_count()
{
  return conv_code_size (ConvBlockType::a, Params::payload_size) * Params::frames_per_bit;
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
  const int frame_count = mark_data_frame_count();
  vector<MixEntry> mix_entries (frame_count * Params::bands_per_frame);

  UpDownGen up_down_gen (Random::Stream::data_up_down);
  int entry = 0;
  for (int f = 0; f < frame_count; f++)
    {
      const int index = data_frame_pos (f);
      UpDownArray up, down;
      up_down_gen.get (f, up, down);

      assert (up.size() == down.size());
      for (size_t i = 0; i < up.size(); i++)
        mix_entries[entry++] = { index, up[i], down[i] };
    }
  Random random (/* seed */ 0, Random::Stream::mix);
  random.shuffle (mix_entries);

  return mix_entries;
}

enum class FrameMod : uint8_t {
  KEEP = 0,
  UP,
  DOWN
};

void
prepare_frame_mod (UpDownGen& up_down_gen, int f, vector<FrameMod>& frame_mod, int data_bit)
{
  UpDownArray up, down;
  up_down_gen.get (f, up, down);
  for (auto u : up)
    frame_mod[u] = data_bit ? FrameMod::UP : FrameMod::DOWN;

  for (auto d : down)
    frame_mod[d] = data_bit ? FrameMod::DOWN : FrameMod::UP;
}

void
apply_frame_mod (const vector<FrameMod>& frame_mod, const vector<complex<float>>& fft_out, vector<complex<float>>& fft_delta_spect)
{
  const float   min_mag = 1e-7;   // avoid computing pow (0.0, -water_delta) which would be inf
  for (size_t i = 0; i < frame_mod.size(); i++)
    {
      if (frame_mod[i] == FrameMod::KEEP)
        continue;

      int data_bit_sign = (frame_mod[i] == FrameMod::UP) ? 1 : -1;
      /*
       * for up bands, we want do use [for a 1 bit]  (pow (mag, 1 - water_delta))
       *
       * this actually increases the amount of energy because mag is less than 1.0
       */
      const float mag = abs (fft_out[i]);
      if (mag > min_mag)
        {
          const float mag_factor = powf (mag, -Params::water_delta * data_bit_sign);

          fft_delta_spect[i] = fft_out[i] * (mag_factor - 1);
        }
    }
}

void
mark_data (vector<vector<FrameMod>>& frame_mod, const vector<int>& bitvec)
{
  assert (bitvec.size() == mark_data_frame_count() / Params::frames_per_bit);
  assert (frame_mod.size() >= mark_data_frame_count());

  const int frame_count = mark_data_frame_count();

  if (Params::mix)
    {
      vector<MixEntry> mix_entries = gen_mix_entries();

      for (int f = 0; f < frame_count; f++)
        {
          for (size_t frame_b = 0; frame_b < Params::bands_per_frame; frame_b++)
            {
              int b = f * Params::bands_per_frame + frame_b;

              const int data_bit = bitvec[f / Params::frames_per_bit];

              const int u = mix_entries[b].up;
              const int d = mix_entries[b].down;
              const int index = mix_entries[b].frame;

              frame_mod[index][u] = data_bit ? FrameMod::UP : FrameMod::DOWN;
              frame_mod[index][d] = data_bit ? FrameMod::DOWN : FrameMod::UP;
            }
        }
    }
  else
    {
      UpDownGen up_down_gen (Random::Stream::data_up_down);
      // sync block always written in linear order (no mix)
      for (int f = 0; f < frame_count; f++)
        {
          size_t index = data_frame_pos (f);

          prepare_frame_mod (up_down_gen, f, frame_mod[index], bitvec[f / Params::frames_per_bit]);
        }
    }
}

size_t
mark_sync_frame_count()
{
  return Params::sync_bits * Params::sync_frames_per_bit;
}

void
mark_sync (vector<vector<FrameMod>>& frame_mod, int ab)
{
  const int frame_count = mark_sync_frame_count();
  assert (frame_mod.size() >= mark_sync_frame_count());

  UpDownGen up_down_gen (Random::Stream::sync_up_down);

  // sync block always written in linear order (no mix)
  for (int f = 0; f < frame_count; f++)
    {
      size_t index = sync_frame_pos (f);
      int    data_bit = (f / Params::sync_frames_per_bit + ab) & 1; /* write 010101 for a block, 101010 for b block */

      prepare_frame_mod (up_down_gen, f, frame_mod[index], data_bit);
    }
}

void
init_pad_mod_vec (vector<vector<FrameMod>>& pad_mod_vec)
{
  UpDownGen up_down_gen (Random::Stream::pad_up_down);

  for (size_t f = 0; f < Params::frames_pad_start; f++)
    {
      vector<FrameMod> mod (Params::max_band + 1);

      prepare_frame_mod (up_down_gen, f, mod, 0);
      pad_mod_vec.push_back (mod);
    }
}

void
init_frame_mod_vec (vector<vector<FrameMod>>& frame_mod_vec, int ab, const vector<int>& bitvec)
{
  frame_mod_vec.resize (mark_sync_frame_count() + mark_data_frame_count());

  for (auto& frame_mod : frame_mod_vec)
    frame_mod.resize (Params::max_band + 1);

  mark_sync (frame_mod_vec, ab);
  mark_data (frame_mod_vec, bitvec);
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

/* synthesizes a watermark stream (overlap add with synthesis window)
 *
 * input:  per-channel fft delta values (always one frame)
 * output: samples
 */
class WatermarkSynth
{
  const int     n_channels = 0;
  vector<float> window;
  vector<float> synth_samples;
  bool          first_frame = true;

  void
  generate_window()
  {
    window.resize (Params::frame_size * 3);
    for (size_t i = 0; i < window.size(); i++)
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
        window[i] = (cos (tri*M_PI+M_PI)+1) * 0.5;
      }
  }
public:
  WatermarkSynth (int n_channels) :
    n_channels (n_channels)
  {
    generate_window();
    synth_samples.resize (window.size() * n_channels);
  }
  vector<float>
  run (const vector<vector<complex<float>>>& fft_delta_spect)
  {
    const size_t synth_frame_sz = Params::frame_size * n_channels;
    /* move frame 1 and frame 2 to frame 0 and frame 1 */
    std::copy (&synth_samples[synth_frame_sz], &synth_samples[synth_frame_sz * 3], &synth_samples[0]);
    /* zero out frame 2 */
    std::fill (&synth_samples[synth_frame_sz * 2], &synth_samples[synth_frame_sz * 3], 0);
    for (int ch = 0; ch < n_channels; ch++)
      {
        /* mix watermark signal to output frame */
        vector<float> fft_delta_out = ifft (fft_delta_spect[ch]);

        for (int dframe = 0; dframe <= 2; dframe++)
          {
            const int wstart = dframe * Params::frame_size;

            int pos = dframe * Params::frame_size * n_channels + ch;
            for (size_t x = 0; x < Params::frame_size; x++)
              {
                synth_samples[pos] += fft_delta_out[x] * window[wstart + x];
                pos += n_channels;
              }
          }
      }
    if (first_frame)
      {
        first_frame = false;
        return {};
      }
    else
      {
        vector<float> out_samples (synth_samples.begin(), synth_samples.begin() + Params::frame_size * n_channels);
        return out_samples;
      }
  }
};

/* generates a watermark signal
 *
 * input:  original signal samples (always for one complete frame)
 * output: watermark signal (to be mixed to the original sample)
 */
class WatermarkGen
{
  enum State { PAD, WATERMARK } state = State::PAD;

  const int                 n_channels = 0;
  int                       frame_number = 0;
  int                       frame_bound = Params::frames_pad_start;
  int                       ab = 0;

  FFTAnalyzer               fft_analyzer;
  WatermarkSynth            wm_synth;

  vector<int>               bitvec_a;
  vector<int>               bitvec_b;
  vector<vector<FrameMod>>  pad_mod_vec;
  vector<vector<FrameMod>>  frame_mod_vec_a;
  vector<vector<FrameMod>>  frame_mod_vec_b;
public:
  WatermarkGen (int n_channels, const vector<int>& bitvec_a, const vector<int>& bitvec_b) :
    n_channels (n_channels),
    fft_analyzer (n_channels),
    wm_synth (n_channels),
    bitvec_a (bitvec_a),
    bitvec_b (bitvec_b)
  {
    init_pad_mod_vec (pad_mod_vec);
    init_frame_mod_vec (frame_mod_vec_a, 0, bitvec_a);
  }
  vector<float>
  run (const vector<float>& samples)
  {
    assert (samples.size() == Params::frame_size * n_channels);

    vector<vector<complex<float>>> fft_out = fft_analyzer.run_fft (samples, 0);

    vector<vector<complex<float>>> fft_delta_spect;
    for (int ch = 0; ch < n_channels; ch++)
      fft_delta_spect.push_back (vector<complex<float>> (fft_out.back().size()));

    if (state == State::PAD)
      {
        for (int ch = 0; ch < n_channels; ch++)
          apply_frame_mod (pad_mod_vec[frame_number], fft_out[ch], fft_delta_spect[ch]);
      }
    else if (state == State::WATERMARK)
      {
        for (int ch = 0; ch < n_channels; ch++)
          apply_frame_mod (ab ? frame_mod_vec_b[frame_number] : frame_mod_vec_a[frame_number], fft_out[ch], fft_delta_spect[ch]);
      }

    frame_number++;
    if (frame_number == frame_bound)
      {
        frame_number = 0;

        if (state == PAD)
          {
            state = WATERMARK;
            frame_bound = mark_sync_frame_count() + mark_data_frame_count();
          }
        else if (state == WATERMARK)
          {
            ab = (ab + 1) & 1; // write A|B|A|B|...
            frame_bound = mark_sync_frame_count() + mark_data_frame_count();

            if (frame_mod_vec_b.empty())
              {
                // we initialize this only when we need it to minimize startup latency
                init_frame_mod_vec (frame_mod_vec_b, 1, bitvec_b);
              }
          }
      }
    return wm_synth.run (fft_delta_spect);
  }
};

class AudioBuffer
{
  const int     n_channels = 0;
  vector<float> buffer;

public:
  AudioBuffer (int n_channels) :
    n_channels (n_channels)
  {
  }
  void
  write_frames (const vector<float>& samples)
  {
    buffer.insert (buffer.end(), samples.begin(), samples.end());
  }
  vector<float>
  read_frames (size_t frames)
  {
    assert (frames * n_channels <= buffer.size());
    const auto begin = buffer.begin();
    const auto end   = begin + frames * n_channels;
    vector<float> result (begin, end);
    buffer.erase (begin, end);
    return result;
  }
  size_t
  can_read_frames() const
  {
    return buffer.size() / n_channels;
  }
};

class ResamplerImpl
{
public:
  virtual
  ~ResamplerImpl()
  {
  }

  virtual void          write_frames (const vector<float>& frames) = 0;
  virtual vector<float> read_frames (size_t frames) = 0;
  virtual size_t        can_read_frames() const = 0;
};

template<class Resampler>
class BufferedResamplerImpl : public ResamplerImpl
{
  const int     n_channels = 0;
  bool          first_write = true;
  Resampler     m_resampler;

  vector<float> buffer;
public:
  BufferedResamplerImpl (int n_channels) :
    n_channels (n_channels)
  {
  }
  Resampler&
  resampler()
  {
    return m_resampler;
  }
  void
  write_frames (const vector<float>& frames)
  {
    if (first_write)
      {
        /* avoid timeshift: zita needs k/2 - 1 samples before the actual input */
        m_resampler.inp_count = m_resampler.inpsize () / 2 - 1;
        m_resampler.inp_data  = nullptr;

        m_resampler.out_count = 1000000; // <- just needs to be large enough that all input is consumed
        m_resampler.out_data  = nullptr;
        m_resampler.process();

        first_write = false;
      }

    uint start = 0;
    do
      {
        const int out_count = Params::frame_size;
        float out[out_count * n_channels];

        m_resampler.out_count = out_count;
        m_resampler.out_data  = out;

        m_resampler.inp_count = frames.size() / n_channels - start;
        m_resampler.inp_data  = const_cast<float *> (&frames[start * n_channels]);
        m_resampler.process();

        size_t count = out_count - m_resampler.out_count;
        buffer.insert (buffer.end(), out, out + count * n_channels);

        start += frames.size() / n_channels - start - m_resampler.inp_count;
      }
    while (start != frames.size() / n_channels);
  }
  vector<float>
  read_frames (size_t frames)
  {
    assert (frames * n_channels <= buffer.size());
    const auto begin = buffer.begin();
    const auto end   = begin + frames * n_channels;
    vector<float> result (begin, end);
    buffer.erase (begin, end);
    return result;
  }
  size_t
  can_read_frames() const
  {
    return buffer.size() / n_channels;
  }
};

ResamplerImpl *
create_resampler (int n_channels, int old_rate, int new_rate)
{
  if (old_rate == new_rate)
    {
      return nullptr; // should not be using create_resampler for that case
    }
  else
    {
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
      const int hlen = 16;

      auto resampler = new BufferedResamplerImpl<Resampler> (n_channels);
      if (resampler->resampler().setup (old_rate, new_rate, n_channels, hlen) == 0)
        {
          return resampler;
        }
      else
        delete resampler;

      auto vresampler = new BufferedResamplerImpl<VResampler> (n_channels);
      const double ratio = double (new_rate) / old_rate;
      if (vresampler->resampler().setup (ratio, n_channels, hlen) == 0)
        {
          return vresampler;
        }
      else
        {
          error ("audiowmark: resampling from old_rate=%d to new_rate=%d not implemented\n", old_rate, new_rate);
          delete vresampler;
          return nullptr;
        }
    }
}

/* generate a watermark at Params::mark_sample_rate and resample to whatever the original signal has
 *
 * input:  samples from original signal (always one frame)
 * output: watermark signal resampled to original signal sample rate
 */
class WatermarkResampler
{
  std::unique_ptr<ResamplerImpl> in_resampler;
  std::unique_ptr<ResamplerImpl> out_resampler;
  WatermarkGen                   wm_gen;
  bool                           need_resampler = false;
public:
  WatermarkResampler (int n_channels, int input_rate, vector<int>& bitvec_a, vector<int>& bitvec_b) :
    wm_gen (n_channels , bitvec_a, bitvec_b)
  {
    need_resampler = (input_rate != Params::mark_sample_rate);

    if (need_resampler)
      {
        in_resampler.reset (create_resampler (n_channels, input_rate, Params::mark_sample_rate));
        out_resampler.reset (create_resampler (n_channels, Params::mark_sample_rate, input_rate));
      }
  }
  bool
  init_ok()
  {
    if (need_resampler)
      return (in_resampler && out_resampler);
    else
      return true;
  }
  vector<float>
  run (const vector<float>& samples)
  {
    if (!need_resampler)
      {
        /* cheap case: if no resampling is necessary, just generate the watermark signal */
        return wm_gen.run (samples);
      }

    /* resample to the watermark sample rate */
    in_resampler->write_frames (samples);
    while (in_resampler->can_read_frames() >= Params::frame_size)
      {
        vector<float> r_samples = in_resampler->read_frames (Params::frame_size);

        /* generate watermark at normalized sample rate */
        vector<float> wm_samples = wm_gen.run (r_samples);

        /* resample back to the original sample rate of the audio file */
        out_resampler->write_frames (wm_samples);
      }

    size_t to_read = out_resampler->can_read_frames();
    return out_resampler->read_frames (to_read);
  }
};

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
  info ("Input:        %s\n", infile.c_str());
  info ("Output:       %s\n", outfile.c_str());
  info ("Message:      %s\n", bit_vec_to_str (bitvec).c_str());
  info ("Strength:     %.6g\n\n", Params::water_delta * 1000);

  /* add forward error correction, bitvec will now be a lot larger */
  auto bitvec_a = randomize_bit_order (conv_encode (ConvBlockType::a, bitvec), /* encode */ true);
  auto bitvec_b = randomize_bit_order (conv_encode (ConvBlockType::b, bitvec), /* encode */ true);

  auto in_stream = std::make_unique<SFInputStream> (); // FIXME: need virtual constructor
  if (!in_stream->open (infile))
    {
      fprintf (stderr, "audiowmark: error opening %s: %s\n", infile.c_str(), in_stream->error_blurb());
      return 1;
    }
  int orig_seconds = in_stream->n_frames() / in_stream->sample_rate();
  info ("Time:         %d:%02d\n", orig_seconds / 60, orig_seconds % 60);
  info ("Sample Rate:  %d\n", in_stream->sample_rate());
  info ("Channels:     %d\n", in_stream->n_channels());

  const int out_bit_depth = in_stream->bit_depth() > 16 ? 24 : 16;
  std::unique_ptr<AudioOutputStream> out_stream;
  if (outfile == "-")
    {
      StdoutWavOutputStream *swstream = new StdoutWavOutputStream();
      out_stream.reset (swstream);
      if (!swstream->open (in_stream->n_channels(), in_stream->sample_rate(), out_bit_depth, in_stream->n_frames()))
        {
          fprintf (stderr, "audiowmark: error writing to -\n"); //%s: %s\n", outfile.c_str(), out_wav_data.error_blurb()); FIXME
          return 1;
        }
    }
  else
    {
      SFOutputStream *sfostream = new SFOutputStream();
      out_stream.reset (sfostream);
      if (!sfostream->open (outfile, in_stream->n_channels(), in_stream->sample_rate(), out_bit_depth, in_stream->n_frames()))
        {
          error ("audiowmark: error writing to %s\n", outfile.c_str()); // FIXME, sfostream->error_blurb());
          return 1;
        }
    }
  vector<float> samples;

  const int n_channels = in_stream->n_channels();
  AudioBuffer audio_buffer (n_channels);
  WatermarkResampler wm_resampler (n_channels, in_stream->sample_rate(), bitvec_a, bitvec_b);
  if (!wm_resampler.init_ok())
    return 1;

  size_t total_input_frames = 0;
  size_t total_output_frames = 0;
  while (true)
    {
      Error err = Error::Code::NONE;

      err = in_stream->read_frames (samples, Params::frame_size);
      if (err)
        {
          error ("audiowmark: input stream read failed: %s\n", err.message());
          return 1;
        }
      total_input_frames += samples.size() / n_channels;

      if (samples.size() < Params::frame_size * n_channels)
        {
          if (total_input_frames == total_output_frames)
            break;

          /* zero sample padding after the actual input */
          samples.resize (Params::frame_size * n_channels);
        }
      audio_buffer.write_frames (samples);
      samples = wm_resampler.run (samples);
      size_t to_read = samples.size() / n_channels;
      vector<float> orig_samples  = audio_buffer.read_frames (to_read);
      assert (samples.size() == orig_samples.size());

      for (size_t i = 0; i < samples.size(); i++)
        samples[i] += orig_samples[i];

      size_t max_write_frames = total_input_frames - total_output_frames;
      if (samples.size() > max_write_frames * n_channels)
        samples.resize (max_write_frames * n_channels);

      err = out_stream->write_frames (samples);
      if (err)
        {
          error ("audiowmark output write failed: %s\n", err.message());
          return 1;
        }
      total_output_frames += samples.size() / n_channels;
    }
#if 0
  if (Params::snr)
    {
      /* compute/show signal to noise ratio */
      double delta_power = 0;
      double signal_power = 0;
      for (size_t i = 0; i < samples.size(); i++)
        {
          const double orig_scaled = samples[i];      // original sample
          const double delta       = out_signal[i];   // watermark

          delta_power += delta * delta;
          signal_power += orig_scaled * orig_scaled;
        }
      delta_power /= samples.size();
      signal_power /= samples.size();

      printf ("SNR:          %f dB\n", 10 * log10 (signal_power / delta_power));
    }
  float max_value = 1e-6;
  for (size_t i = 0; i < samples.size(); i++)
    {
      /* Typically the original samples are already in range [-1;1]. However in
       * some cases (mp3 loader), the samples are not fully normalized; in those
       * cases, for volume normalization we treat them as-if they had been
       * clipped already; final clipping will be done while saving.
       */
      const float x = bound<float> (-1, samples[i], 1);
      const float value = fabsf (x + out_signal[i]);
      if (value > max_value)
        max_value = value;
    }

  // scale (samples + watermark) down if necessary to avoid clipping
  const float scale = min (1.0 / max_value, 1.0);
  for (size_t i = 0; i < samples.size(); i++)
    samples[i] = (samples[i] + out_signal[i]) * scale;

  printf ("Data Blocks:  %d\n", data_blocks);
  printf ("Volume Norm:  %.3f (%.2f dB)\n", scale, db_from_factor (scale, -96));
#endif
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
mix_decode (vector<vector<complex<float>>>& fft_out, int n_channels)
{
  vector<float> raw_bit_vec;

  const int frame_count = mark_data_frame_count();

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
linear_decode (vector<vector<complex<float>>>& fft_out, int n_channels)
{
  UpDownGen     up_down_gen (Random::Stream::data_up_down);
  vector<float> raw_bit_vec;

  const int frame_count = mark_data_frame_count();

  double umag = 0, dmag = 0;
  for (int f = 0; f < frame_count; f++)
    {
      for (int ch = 0; ch < n_channels; ch++)
        {
          const size_t index = data_frame_pos (f) * n_channels + ch;
          UpDownArray up, down;
          up_down_gen.get (f, up, down);

          const double min_db = -96;
          for (auto u : up)
            umag += db_from_factor (abs (fft_out[index][u]), min_db);

          for (auto d : down)
            dmag += db_from_factor (abs (fft_out[index][d]), min_db);
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

    UpDownGen up_down_gen (Random::Stream::sync_up_down);
    size_t n_bands = Params::max_band - Params::min_band + 1;
    for (int bit = 0; bit < Params::sync_bits; bit++)
      {
        for (int f = 0; f < Params::sync_frames_per_bit; f++)
          {
            UpDownArray frame_up, frame_down;
            up_down_gen.get (f + bit * Params::sync_frames_per_bit, frame_up, frame_down);

            for (auto u : frame_up)
              up[bit].push_back (u - Params::min_band + sync_frame_pos (f + bit * Params::sync_frames_per_bit) * n_bands * wav_data.n_channels());

            for (auto d : frame_down)
              down[bit].push_back (d - Params::min_band + sync_frame_pos (f + bit * Params::sync_frames_per_bit) * n_bands * wav_data.n_channels());
          }
        sort (up[bit].begin(), up[bit].end());
        sort (down[bit].begin(), down[bit].end());
      }
  }
  double
  sync_decode (const WavData& wav_data, const size_t start_frame, const vector<float>& fft_out_db, ConvBlockType *block_type)
  {
    double sync_quality = 0;

    size_t n_bands = Params::max_band - Params::min_band + 1;
    for (int bit = 0; bit < Params::sync_bits; bit++)
      {
        float umag = 0, dmag = 0;

        for (int ch = 0; ch < wav_data.n_channels(); ch++)
          {
            const int index = (start_frame * wav_data.n_channels() + ch) * n_bands;

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

    if (sync_quality < 0)
      {
        *block_type = ConvBlockType::b;
        return -sync_quality;
      }
    else
      {
        *block_type = ConvBlockType::a;
        return sync_quality;
      }
  }
public:
  struct Score {
    size_t        index;
    double        quality;
    ConvBlockType block_type;
  };
  vector<Score>
  search (const WavData& wav_data)
  {
    vector<Score> result_scores;
    vector<Score> sync_scores;

    if (Params::test_no_sync)
      {
        const size_t expect0 = Params::frames_pad_start * Params::frame_size;
        const size_t expect_step = (mark_sync_frame_count() + mark_data_frame_count()) * Params::frame_size;
        const size_t expect_end = frame_count (wav_data) * Params::frame_size;

        int ab = 0;
        for (size_t expect_index = expect0; expect_index + expect_step < expect_end; expect_index += expect_step)
          result_scores.push_back (Score { expect_index, 1.0, (ab++ & 1) ? ConvBlockType::b : ConvBlockType::a });

        return result_scores;
      }
    init_up_down (wav_data);

    vector<float> fft_db;

    // compute multiple time-shifted fft vectors
    size_t n_bands = Params::max_band - Params::min_band + 1;
    for (size_t sync_shift = 0; sync_shift < Params::frame_size; sync_shift += Params::sync_search_step)
      {
        sync_fft (wav_data, sync_shift, frame_count (wav_data) - 1, fft_db, /* want all frames */ {});
        for (int start_frame = 0; start_frame < frame_count (wav_data); start_frame++)
          {
            const size_t sync_index = start_frame * Params::frame_size + sync_shift;
            if ((start_frame + mark_sync_frame_count() + mark_data_frame_count()) * wav_data.n_channels() * n_bands < fft_db.size())
              {
                ConvBlockType block_type;
                double quality = sync_decode (wav_data, start_frame, fft_db, &block_type);
                // printf ("%zd %f\n", sync_index, quality);
                sync_scores.emplace_back (Score { sync_index, quality, block_type });
              }
          }
      }
    sort (sync_scores.begin(), sync_scores.end(), [] (const Score& a, const Score &b) { return a.index < b.index; });

    vector<int> want_frames (mark_sync_frame_count() + mark_data_frame_count());
    for (size_t f = 0; f < mark_sync_frame_count(); f++)
      want_frames[sync_frame_pos (f)] = 1;

    /* for strength 8 and above:
     *   -> more false positive candidates are rejected, so we can use a lower threshold
     *
     * for strength 7 and below:
     *   -> we need a higher threshold, because otherwise watermark detection takes too long
     */
    const double strength = Params::water_delta * 1000;
    const double sync_threshold1 = strength > 7.5 ? 0.4 : 0.5;

    for (size_t i = 0; i < sync_scores.size(); i++)
      {
        // printf ("%zd %f\n", sync_scores[i].index, sync_scores[i].quality);
        if (sync_scores[i].quality > sync_threshold1)
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
                double best_quality       = sync_scores[i].quality;
                size_t best_index         = sync_scores[i].index;
                ConvBlockType best_block_type = sync_scores[i].block_type; /* doesn't really change during refinement */

                int start = std::max (int (sync_scores[i].index) - Params::sync_search_step, 0);
                int end   = sync_scores[i].index + Params::sync_search_step;
                for (int fine_index = start; fine_index <= end; fine_index += Params::sync_search_fine)
                  {
                    sync_fft (wav_data, fine_index, mark_sync_frame_count() + mark_data_frame_count(), fft_db, want_frames);
                    if (fft_db.size())
                      {
                        ConvBlockType block_type;
                        double        q = sync_decode (wav_data, 0, fft_db, &block_type);

                        if (q > best_quality)
                          {
                            best_quality = q;
                            best_index   = fine_index;
                          }
                      }
                  }
                //printf (" => refined: %zd %s %f\n", best_index, find_closest_sync (best_index), best_quality);
                if (best_quality > Params::sync_threshold2)
                  result_scores.push_back (Score { best_index, best_quality, best_block_type });
              }
          }
      }
    return result_scores;
  }
private:
  void
  sync_fft (const WavData& wav_data, size_t index, size_t frame_count, vector<float>& fft_out_db, const vector<int>& want_frames)
  {
    fft_out_db.clear();

    /* read past end? -> fail */
    if (wav_data.n_values() < (index + frame_count * Params::frame_size) * wav_data.n_channels())
      return;

    FFTAnalyzer fft_analyzer (wav_data.n_channels());
    const vector<float>& samples = wav_data.samples();
    const size_t n_bands = Params::max_band - Params::min_band + 1;
    int out_pos = 0;

    fft_out_db.resize (wav_data.n_channels() * n_bands * frame_count);

    for (size_t f = 0; f < frame_count; f++)
      {
        const double min_db = -96;
        if (want_frames.size() && !want_frames[f])
          {
            for (int ch = 0; ch < wav_data.n_channels(); ch++)
              for (int i = Params::min_band; i <= Params::max_band; i++)
                fft_out_db[out_pos++] = min_db;
          }
        else
          {
            vector<vector<complex<float>>> frame_result = fft_analyzer.run_fft (samples, index + f * Params::frame_size);

            /* computing db-magnitude is expensive, so we better do it here */
            for (int ch = 0; ch < wav_data.n_channels(); ch++)
              for (int i = Params::min_band; i <= Params::max_band; i++)
                fft_out_db[out_pos++] = db_from_factor (abs (frame_result[ch][i]), min_db);
          }
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

  auto report_pattern = [&] (SyncFinder::Score sync_score, const vector<int>& bit_vec, float decode_error)
  {
    if (sync_score.index)
      {
        const char *block_str = nullptr;

        switch (sync_score.block_type)
          {
            case ConvBlockType::a:  block_str = "A";
                                    break;
            case ConvBlockType::b:  block_str = "B";
                                    break;
            case ConvBlockType::ab: block_str = "AB";
                                    break;
          }
        const int seconds = sync_score.index / wav_data.sample_rate();
        printf ("pattern %2d:%02d %s %.3f %.3f %s\n", seconds / 60, seconds % 60, bit_vec_to_str (bit_vec).c_str(),
                sync_score.quality, decode_error, block_str);
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

  vector<float> raw_bit_vec_all (conv_code_size (ConvBlockType::ab, Params::payload_size));
  vector<int>   raw_bit_vec_norm (2);

  SyncFinder::Score score_all { 0, 0 };
  SyncFinder::Score score_ab  { 0, 0, ConvBlockType::ab };

  ConvBlockType last_block_type = ConvBlockType::b;
  vector<vector<float>> ab_raw_bit_vec (2);
  vector<float>         ab_quality (2);
  FFTAnalyzer           fft_analyzer (wav_data.n_channels());
  for (auto sync_score : sync_scores)
    {
      const size_t count = mark_sync_frame_count() + mark_data_frame_count();
      const size_t index = sync_score.index;
      const int    ab = (sync_score.block_type == ConvBlockType::b); /* A -> 0, B -> 1 */

      auto fft_range_out = fft_analyzer.fft_range (wav_data.samples(), index, count);
      if (fft_range_out.size())
        {
          /* ---- retrieve bits from watermark ---- */
          vector<float> raw_bit_vec;
          if (Params::mix)
            {
              raw_bit_vec = mix_decode (fft_range_out, wav_data.n_channels());
            }
          else
            {
              raw_bit_vec = linear_decode (fft_range_out, wav_data.n_channels());
            }
          assert (raw_bit_vec.size() == conv_code_size (ConvBlockType::a, Params::payload_size));

          raw_bit_vec = randomize_bit_order (raw_bit_vec, /* encode */ false);

          /* ---- deal with this pattern ---- */
          float decode_error = 0;
          vector<int> bit_vec = conv_decode_soft (sync_score.block_type, normalize_soft_bits (raw_bit_vec), &decode_error);

          report_pattern (sync_score, bit_vec, decode_error);

          /* ---- update "all" pattern ---- */
          score_all.quality += sync_score.quality;

          for (size_t i = 0; i < raw_bit_vec.size(); i++)
            {
              raw_bit_vec_all[i * 2 + ab] += raw_bit_vec[i];
            }
          raw_bit_vec_norm[ab]++;

          /* ---- if last block was A & this block is B => deal with combined AB block */
          ab_raw_bit_vec[ab] = raw_bit_vec;
          ab_quality[ab]     = sync_score.quality;
          if (last_block_type == ConvBlockType::a && sync_score.block_type == ConvBlockType::b)
            {
              /* join A and B block -> AB block */
              vector<float> ab_bits (raw_bit_vec.size() * 2);
              for (size_t i = 0; i <  raw_bit_vec.size(); i++)
                {
                  ab_bits[i * 2] = ab_raw_bit_vec[0][i];
                  ab_bits[i * 2 + 1] = ab_raw_bit_vec[1][i];
                }
              vector<int> bit_vec = conv_decode_soft (ConvBlockType::ab, normalize_soft_bits (ab_bits), &decode_error);
              score_ab.index = sync_score.index;
              score_ab.quality = (ab_quality[0] + ab_quality[1]) / 2;
              report_pattern (score_ab, bit_vec, decode_error);
            }
          last_block_type = sync_score.block_type;
        }
    }
  if (total_count > 1) /* all pattern: average soft bits of all watermarks and decode */
    {
      for (size_t i = 0; i < raw_bit_vec_all.size(); i += 2)
        {
          raw_bit_vec_all[i]     /= max (raw_bit_vec_norm[0], 1); /* normalize A soft bits with number of A blocks */
          raw_bit_vec_all[i + 1] /= max (raw_bit_vec_norm[1], 1); /* normalize B soft bits with number of B blocks */
        }
      score_all.quality /= raw_bit_vec_norm[0] + raw_bit_vec_norm[1];

      vector<float> soft_bit_vec = normalize_soft_bits (raw_bit_vec_all);

      float decode_error = 0;
      vector<int> bit_vec = conv_decode_soft (ConvBlockType::ab, soft_bit_vec, &decode_error);

      report_pattern (score_all, bit_vec, decode_error);
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

  if (Params::test_truncate)
    {
      const size_t  want_n_samples = wav_data.sample_rate() * wav_data.n_channels() * Params::test_truncate;
      vector<float> short_samples  = wav_data.samples();

      if (want_n_samples < short_samples.size())
        {
          short_samples.resize (want_n_samples);
          wav_data.set_samples (short_samples);
        }
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
test_subtract (const string& infile1, const string& infile2, const string& outfile)
{
  WavData in1_data;
  if (!in1_data.load (infile1))
    {
      error ("audiowmark: error loading %s: %s\n", infile1.c_str(), in1_data.error_blurb());
      return 1;
    }
  WavData in2_data;
  if (!in2_data.load (infile2))
    {
      error ("audiowmark: error loading %s: %s\n", infile2.c_str(), in2_data.error_blurb());
      return 1;
    }
  if (in1_data.n_values() != in2_data.n_values())
    {
      int64_t l1 = in1_data.n_values();
      int64_t l2 = in2_data.n_values();
      size_t  delta = std::abs (l1 - l2);
      warning ("audiowmark: size mismatch: %zd frames\n", delta / in1_data.n_channels());
      warning (" - %s frames: %zd\n", infile1.c_str(), in1_data.n_values() / in1_data.n_channels());
      warning (" - %s frames: %zd\n", infile2.c_str(), in2_data.n_values() / in2_data.n_channels());
    }
  assert (in1_data.n_channels() == in2_data.n_channels());

  const auto& in1_signal = in1_data.samples();
  const auto& in2_signal = in2_data.samples();
  size_t len = std::min (in1_data.n_values(), in2_data.n_values());
  vector<float> out_signal;
  for (size_t i = 0; i < len; i++)
    out_signal.push_back (in1_signal[i] - in2_signal[i]);

  WavData out_wav_data (out_signal, in1_data.n_channels(), in1_data.sample_rate(), in1_data.bit_depth());
  if (!out_wav_data.save (outfile))
    {
      fprintf (stderr, "audiowmark: error saving %s: %s\n", outfile.c_str(), out_wav_data.error_blurb());
      return 1;
    }
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
  else if (op == "cut-start" && argc == 5)
    {
      cut_start (argv[2], argv[3], argv[4]);
    }
  else if (op == "test-subtract" && argc == 5)
    {
      test_subtract (argv[2], argv[3], argv[4]);
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
