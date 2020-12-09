/*
 * Copyright (C) 2018-2020 Stefan Westerfeld
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "wmspeed.hh"
#include "wmcommon.hh"
#include "syncfinder.hh"
#include "threadpool.hh"
#include "fft.hh"

#include <algorithm>
#include <zita-resampler/vresampler.h>

using std::vector;
using std::sort;

/* FIXME: dedup this template */
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

static WavData
resample_ratio (const WavData& wav_data, double ratio, int new_rate)
{
  const int hlen = 16;
  const vector<float>& in = wav_data.samples();
  vector<float> out (lrint (in.size() / wav_data.n_channels() * ratio) * wav_data.n_channels());

  VResampler vresampler;
  if (vresampler.setup (ratio, wav_data.n_channels(), hlen) != 0)
    {
      error ("audiowmark: failed to setup vresampler with ratio=%f\n", ratio);
      exit (1);
    }

  process_resampler (vresampler, in, out);
  return WavData (out, wav_data.n_channels(), new_rate, wav_data.bit_depth());
}

static WavData
truncate (const WavData& in_data, double seconds)
{
  WavData out_data = in_data;
  vector<float> short_samples = in_data.samples();

  const size_t want_n_samples = lrint (in_data.sample_rate() * seconds) * in_data.n_channels();
  if (short_samples.size() > want_n_samples)
    {
      short_samples.resize (want_n_samples);
      out_data.set_samples (short_samples);
    }
  return out_data;
}

class SpeedSync
{
public:
  struct Score
  {
    double speed = 0;
    double quality = 0;
  };
private:
  struct Mags {
    float umag = 0;
    float dmag = 0;
  };
  vector<vector<SyncFinder::FrameBit>> sync_bits;
  vector<vector<Mags>> fft_sync_bits;

  void prepare_mags();
  void compare (double relative_speed);

  std::mutex mutex;
  vector<Score> result_scores;
  const WavData& in_data;
  const double center;
  const double step;
  const int    n_steps;
  const double seconds;
  const int    frames_per_block;
public:
  SpeedSync (const WavData& in_data, double center, double step, int n_steps, double seconds) :
    in_data (in_data),
    center (center),
    step (step),
    n_steps (n_steps),
    seconds (seconds),
    frames_per_block (mark_sync_frame_count() + mark_data_frame_count())
  {
    // constructor is run in the main thread; everything that is not thread-safe must happen here
    SyncFinder sync_finder;

    sync_bits = sync_finder.get_sync_bits (in_data, SyncFinder::Mode::BLOCK);
  }
  void
  prepare_job (ThreadPool& thread_pool)
  {
    thread_pool.add_job ([this]() { prepare_mags(); });
  }

  void
  search (ThreadPool& thread_pool)
  {
    for (int p = -n_steps; p <= n_steps; p++)
      {
        const double relative_speed = pow (step, p);

        thread_pool.add_job ([relative_speed, this]() { compare (relative_speed); });
      }
  }

  vector<Score>
  get_scores()
  {
    return result_scores;
  }
};

struct SpeedScanParams
{
  double seconds        = 0;
  double step           = 0;
  int    n_steps        = 0;
  int    n_center_steps = 0;
};

static double
speed_scan (ThreadPool& thread_pool, const WavData& in_data, const SpeedScanParams& params, double speed)
{
  vector<SpeedSync::Score> scores;

  /* n_center_steps / n_steps settings: speed approximately 0.8..1.25 */

  vector<std::unique_ptr<SpeedSync>> speed_sync;

  auto t = get_time();
  for (int c = -params.n_center_steps; c <= params.n_center_steps; c++)
    {
      double c_speed = speed * pow (params.step, c * (params.n_steps * 2 + 1));

      speed_sync.push_back (std::make_unique<SpeedSync> (in_data, c_speed, params.step, params.n_steps, params.seconds));
    }

  for (auto& s : speed_sync)
    s->prepare_job (thread_pool);

  thread_pool.wait_all();
  printf ("## wait prepare jobs: %f\n", get_time() - t);
  t=get_time();

  for (auto& s : speed_sync)
    s->search (thread_pool);

  thread_pool.wait_all();
  printf ("## wait search jobs: %f\n", get_time() - t);
  t=get_time();

  for (auto& s : speed_sync)
    {
      vector<SpeedSync::Score> step_scores = s->get_scores();
      scores.insert (scores.end(), step_scores.begin(), step_scores.end());
    }
  sort (scores.begin(), scores.end(), [] (SpeedSync::Score s_a, SpeedSync::Score s_b)
    {
      if (s_a.quality == s_b.quality)
        return s_a.speed > s_b.speed;
       return s_a.quality > s_b.quality;
    });

  // we could search the N best matches, but using the best result works well in practice
  SpeedSync::Score best_s = scores[0];

  printf ("detect_speed_%.0f %f %f %f\n",
    params.seconds,
    best_s.speed,
    best_s.quality,
    100 * fabs (best_s.speed - Params::detect_speed_hint) / Params::detect_speed_hint);
  return best_s.speed;
}

/* FIXME: is this the best choice */
inline double
window_hamming (double x) /* sharp (rectangle) cutoffs at boundaries */
{
  if (fabs (x) > 1)
    return 0;

  return 0.54 + 0.46 * cos (M_PI * x);
}

void
SpeedSync::prepare_mags()
{
  WavData in_data_trc (truncate (in_data, seconds / center));

  // we downsample the audio by factor 2 to improve performance
  WavData in_data_sub (resample_ratio (in_data_trc, center / 2, Params::mark_sample_rate / 2));

  const int sub_frame_size = Params::frame_size / 2;
  const int sub_sync_search_step = Params::sync_search_step / 2;

  double window_weight = 0;
  float window[sub_frame_size];
  for (size_t i = 0; i < sub_frame_size; i++)
    {
      const double fsize_2 = sub_frame_size / 2.0;
      // const double win =  window_cos ((i - fsize_2) / fsize_2);
      const double win = window_hamming ((i - fsize_2) / fsize_2);
      //const double win = 1;
      window[i] = win;
      window_weight += win;
    }

  /* normalize window using window weight */
  for (size_t i = 0; i < sub_frame_size; i++)
    {
      window[i] *= 2.0 / window_weight;
    }

  FFTProcessor fft_processor (sub_frame_size);

  float *in = fft_processor.in();
  float *out = fft_processor.out();

  fft_sync_bits.clear();
  size_t pos = 0;
  while (pos + sub_frame_size < in_data_sub.n_frames())
    {
      const std::vector<float>& samples = in_data_sub.samples();
      vector<float> fft_out_db;

      for (int ch = 0; ch < in_data_sub.n_channels(); ch++)
        {
          for (int i = 0; i < sub_frame_size; i++)
            {
              in[i] = samples[ch + (pos + i) * in_data_sub.n_channels()] * window[i];
            }
          fft_processor.fft();

          for (int i = Params::min_band; i <= Params::max_band; i++)
            {
              const float min_db = -96;
              float re = out[i * 2];
              float im = out[i * 2 + 1];
              float abs = sqrt (re * re + im * im); // FIXME: could avoid sqrt here
              fft_out_db.push_back (db_from_factor (abs, min_db));
            }
        }
      vector<Mags> mags;
      for (size_t bit = 0; bit < sync_bits.size(); bit++)
        {
          const vector<SyncFinder::FrameBit>& frame_bits = sync_bits[bit];
          for (const auto& frame_bit : frame_bits)
            {
              float umag = 0, dmag = 0;

              for (size_t i = 0; i < frame_bit.up.size(); i++)
                {
                  umag += fft_out_db[frame_bit.up[i]];
                  dmag += fft_out_db[frame_bit.down[i]];
                }
              mags.push_back (Mags {umag, dmag});
            }
        }
      fft_sync_bits.push_back (mags);
      pos += sub_sync_search_step;
    }
}

void
SpeedSync::compare (double relative_speed)
{
  const int pad_start = frames_per_block * /* HACK */ 4;
  Score best_score;
  // FIXME: pad_start must be scaled with speed
  for (int offset = -pad_start; offset < 0; offset++)
    {
      double sync_quality = 0;
      int mi = 0;
      int frame_bit_count = 0;
      int bit_count = 0;
      for (size_t sync_bit = 0; sync_bit < Params::sync_bits; sync_bit++)
        {
          const vector<SyncFinder::FrameBit>& frame_bits = sync_bits[sync_bit];

          float umag = 0;
          float dmag = 0;
          for (size_t f = 0; f < Params::sync_frames_per_bit; f++)
            {
              const int index1 = (offset + frame_bits[f].frame * /* HACK */ 4) / relative_speed;
              if (index1 >= 0 && index1 < (int) fft_sync_bits.size())
                {
                  umag += fft_sync_bits[index1][mi].umag;
                  dmag += fft_sync_bits[index1][mi].dmag;
                  frame_bit_count++;
                }
              // FIXME: probably better compute double index
              const int index2 = index1 + (frames_per_block * /* HACK */ 4) / relative_speed;
              if (index2 >= 0 && index2 < (int) fft_sync_bits.size())
                {
                  umag += fft_sync_bits[index2][mi].dmag;
                  dmag += fft_sync_bits[index2][mi].umag;
                  frame_bit_count++;
                }
              mi++;
            }
          /* convert avoiding bias, raw_bit < 0 => 0 bit received; raw_bit > 0 => 1 bit received */
          double raw_bit;
          if (umag == 0 || dmag == 0)
            {
              raw_bit = 0;
            }
          else if (umag < dmag)
            {
              raw_bit = 1 - umag / dmag;
            }
          else
            {
              raw_bit = dmag / umag - 1;
            }
          const int expect_data_bit = sync_bit & 1; /* expect 010101 */
          const double q = expect_data_bit ? raw_bit : -raw_bit;
          sync_quality += q * frame_bit_count;
          bit_count += frame_bit_count;
        }
      if (bit_count)
        {
          sync_quality /= bit_count;
          sync_quality = fabs (sync_quality);
          //sync_quality = fabs (sync_finder.normalize_sync_quality (sync_quality));
          //printf ("%d %f\n", offset, fabs (sync_quality));
          if (sync_quality > best_score.quality)
            {
              best_score.quality = sync_quality;
              best_score.speed = relative_speed * center;
            }
        }
    }
  //printf ("%f %f\n", best_score.speed, best_score.quality);
  std::lock_guard<std::mutex> lg (mutex);
  result_scores.push_back (best_score);
}

static double
get_clip_location (const WavData& in_data)
{
  Random rng (0, Random::Stream::speed_clip);

  /* to improve performance, we don't hash all samples but just a few */
  const vector<float>& samples = in_data.samples();
  vector<float> xsamples;
  for (size_t p = 0; p < samples.size(); p += rng() % 1000)
    xsamples.push_back (samples[p]);

  rng.seed (Random::seed_from_hash (xsamples), Random::Stream::speed_clip);

  return rng.random_double();
}

static WavData
get_speed_clip (double location, const WavData& in_data, double clip_seconds)
{
  double end_sec = double (in_data.n_frames()) / in_data.sample_rate();
  double start_sec = location * (end_sec - clip_seconds);
  if (start_sec < 0)
    start_sec = 0;

  size_t start_point = start_sec * in_data.sample_rate();
  size_t end_point = std::min<size_t> (start_point + clip_seconds * in_data.sample_rate(), in_data.n_frames());
#if 0
  printf ("[%f %f] l%f\n", double (start_point) / in_data.sample_rate(), double (end_point) / in_data.sample_rate(),
                           double (end_point - start_point) / in_data.sample_rate());
#endif
  vector<float> out_signal (in_data.samples().begin() + start_point * in_data.n_channels(),
                            in_data.samples().begin() + end_point * in_data.n_channels());
  WavData clip_data (out_signal, in_data.n_channels(), in_data.sample_rate(), in_data.bit_depth());
  return clip_data;
}

double
detect_speed (const WavData& in_data)
{
  double clip_location = get_clip_location (in_data);

  /* speed is between 0.8 and 1.25, so we use a clip seconds factor of 1.3 to provide enough samples */
  WavData in_clip_short = get_speed_clip (clip_location, in_data, 21 * 1.3);
  WavData in_clip_long  = get_speed_clip (clip_location, in_data, 50 * 1.3);

  ThreadPool thread_pool;

  /* first pass:  find approximation for speed */
  const SpeedScanParams scan1
    {
      .seconds        = 21,
      .step           = 1.0007,
      .n_steps        = 5,
      .n_center_steps = 28
    };
  double speed = speed_scan (thread_pool, in_clip_short, scan1, /* start speed */ 1.0);

  /* second pass: fast refine (not always perfect) */
  const SpeedScanParams scan2
    {
      .seconds        = 50,
      .step           = 1.00005,
      .n_steps        = 20,
      .n_center_steps = 0
    };
  speed = speed_scan (thread_pool, in_clip_long, scan2, speed);
  return speed;
}
