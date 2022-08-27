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
#include "resample.hh"

#include <algorithm>

using std::vector;
using std::sort;
using std::max;

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

struct SpeedScanParams
{
  double seconds        = 0;
  double step           = 0;
  int    n_steps        = 0;
  int    n_center_steps = 0;
};

class MagMatrix
{
public:
  struct Mags
  {
    float umag = 0;
    float dmag = 0;
  };
private:
  vector<Mags> m_data;
  int m_cols = 0;
  int m_rows = 0;
public:
  Mags&
  operator() (int row, int col)
  {
    return m_data[col * m_rows + row];
  }
  void
  resize (int rows, int cols)
  {
    m_rows = rows;
    m_cols = cols;

    /* - don't preserve contents on resize
     * - free unused memory on resize
     */
    vector<Mags> new_data (m_rows * m_cols);
    m_data.swap (new_data);
  }
  int
  rows()
  {
    return m_rows;
  }
};

class SpeedSync
{
public:
  struct Score
  {
    double speed = 0;
    double quality = 0;
  };
  struct SyncBit
  {
    int bit = 0;
    int frame = 0;
    std::vector<int> up;
    std::vector<int> down;
  };
  struct BitValue
  {
    float umag = 0;
    float dmag = 0;
    int count = 0;
  };
  struct CmpState
  {
    int offset = 0;
    BitValue bit_values[Params::sync_bits];
  };
private:
  static constexpr int OFFSET_SHIFT = 16;

  vector<SyncBit> sync_bits;
  MagMatrix sync_matrix;

  void prepare_mags (const SpeedScanParams& scan_params);
  void compare (double relative_speed);
  template<int BLOCK>
  void compare_bits (vector<CmpState>& cmp_states, double relative_speed);

  std::mutex mutex;
  vector<Score> result_scores;
  const WavData& in_data;
  const double center;
  const int    frames_per_block;

public:
  SpeedSync (const WavData& in_data, double center) :
    in_data (in_data),
    center (center),
    frames_per_block (mark_sync_frame_count() + mark_data_frame_count())
  {
    // constructor is run in the main thread; everything that is not thread-safe must happen here
    SyncFinder sync_finder;

    auto sync_finder_bits = sync_finder.get_sync_bits (in_data, SyncFinder::Mode::BLOCK);
    for (size_t bit = 0; bit < sync_finder_bits.size(); bit++)
      {
        for (const auto& frame_bit : sync_finder_bits[bit])
          {
            SyncBit sb;

            sb.bit    = bit;
            sb.frame  = frame_bit.frame;
            sb.up     = frame_bit.up;
            sb.down   = frame_bit.down;

            sync_bits.push_back (sb);
          }
      }
    std::sort (sync_bits.begin(), sync_bits.end(), [](const auto& s1, const auto& s2) { return s1.frame < s2.frame; });
  }
  void
  start_prepare_job (ThreadPool& thread_pool, const SpeedScanParams& scan_params)
  {
    thread_pool.add_job ([this, &scan_params]() { prepare_mags (scan_params); });
  }

  void
  start_search_jobs (ThreadPool& thread_pool, const SpeedScanParams& scan_params, double speed)
  {
    result_scores.clear();

    for (int p = -scan_params.n_steps; p <= scan_params.n_steps; p++)
      {
        const double relative_speed = pow (scan_params.step, p) * speed / center;

        thread_pool.add_job ([relative_speed, this]() { compare (relative_speed); });
      }
  }

  vector<Score>
  get_scores()
  {
    return result_scores;
  }
  double
  center_speed() const
  {
    return center;
  }
};

void
SpeedSync::prepare_mags (const SpeedScanParams& scan_params)
{
  WavData in_data_trc (truncate (in_data, scan_params.seconds / center));

  // we downsample the audio by factor 2 to improve performance
  WavData in_data_sub (resample_ratio (in_data_trc, center / 2, Params::mark_sample_rate / 2));

  const int sub_frame_size = Params::frame_size / 2;
  const int sub_sync_search_step = Params::sync_search_step / 2;

  vector<float> window = FFTAnalyzer::gen_normalized_window (sub_frame_size);

  FFTProcessor fft_processor (sub_frame_size);

  float *in = fft_processor.in();
  float *out = fft_processor.out();

  /* set mag matrix size */
  int n_sync_rows = 0;
  int n_sync_cols = sync_bits.size();
  for (size_t ppos = 0; ppos + sub_frame_size < in_data_sub.n_frames(); ppos += sub_sync_search_step)
    n_sync_rows++;
  sync_matrix.resize (n_sync_rows, n_sync_cols);

  size_t pos = 0;
  int row = 0;
  while (pos + sub_frame_size < in_data_sub.n_frames())
    {
      int col = 0;
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

              fft_out_db.push_back (db_from_complex (out[i * 2], out[i * 2 + 1], min_db));
            }
        }
      for (const auto& sync_bit : sync_bits)
        {
          float umag = 0, dmag = 0;

          for (size_t i = 0; i < sync_bit.up.size(); i++)
            {
              umag += fft_out_db[sync_bit.up[i]];
              dmag += fft_out_db[sync_bit.down[i]];
            }
          sync_matrix (row, col++) = MagMatrix::Mags {umag, dmag};
        }
      assert (col == n_sync_cols);
      row++;
      pos += sub_sync_search_step;
    }
  assert (row == n_sync_rows);
}

template<int BLOCK> void
SpeedSync::compare_bits (vector<CmpState>& cmp_states, double relative_speed)
{
  const int steps_per_frame = Params::frame_size / Params::sync_search_step;
  const double relative_speed_inv = 1 / relative_speed;

  auto begin = cmp_states.end();
  auto end = cmp_states.end();
  for (size_t mi = 0; mi < sync_bits.size(); mi++)
    {
      const int frame_offset = ((BLOCK * frames_per_block + sync_bits[mi].frame) * steps_per_frame * relative_speed_inv + 0.5) * (1 << OFFSET_SHIFT);

      while (begin > cmp_states.begin())
        {
          auto prev = begin - 1;

          /*
           * don't use OFFSET_SHIFT here; just ensure that index is positive
           * to ensure that shifted value will properly round to nearest frame
           * later on
           */
          int index = prev->offset + frame_offset;
          if (index < 0)
            break;

          begin = prev;
        }
      while (end > cmp_states.begin())
        {
          auto prev = end - 1;

          int index = (prev->offset + frame_offset) >> OFFSET_SHIFT;
          if (index < sync_matrix.rows())
            break;

          end = prev;
        }

      for (auto it = begin; it != end; it++)
        {
          int index = (it->offset + frame_offset) >> OFFSET_SHIFT;

          auto& bv = it->bit_values[sync_bits[mi].bit];
          auto mags = sync_matrix (index, mi);
          if (BLOCK & 1)
            {
              bv.umag += mags.dmag;
              bv.dmag += mags.umag;
            }
          else
            {
              bv.umag += mags.umag;
              bv.dmag += mags.dmag;
            }
          bv.count++;
        }
    }
}

void
SpeedSync::compare (double relative_speed)
{
  const int steps_per_frame = Params::frame_size / Params::sync_search_step;
  const int pad_start = frames_per_block * steps_per_frame + /* add a bit of overlap to handle boundaries */ steps_per_frame;

  assert (steps_per_frame * Params::sync_search_step == Params::frame_size);

  vector<CmpState> cmp_states;
  for (int offset =  -pad_start; offset < 0; offset++)
    {
      CmpState cs;
      cs.offset = offset * ((1 << OFFSET_SHIFT) / relative_speed);
      cmp_states.push_back (cs);
    }

  /*
   * we need to compare 3 blocks here:
   *  - one block is necessary because we need to test all offsets (-pad_start..0)
   *  - two more blocks are necessary since speed detection ScanParams uses 50 seconds at most,
   *    and short payload (12 bits) has a block length of slightly over 30 seconds
   */
  compare_bits<0> (cmp_states, relative_speed);
  compare_bits<1> (cmp_states, relative_speed);
  compare_bits<2> (cmp_states, relative_speed);

  Score best_score;
  for (const auto& cs : cmp_states)
    {
      double sync_quality = 0;
      int bit_count = 0;

      for (size_t bit = 0; bit < Params::sync_bits; bit++)
        {
          const auto& bv = cs.bit_values[bit];

          sync_quality += SyncFinder::bit_quality (bv.umag, bv.dmag, bit) * bv.count;
          bit_count += bv.count;
        }
      if (bit_count)
        {
          sync_quality /= bit_count;
          sync_quality = fabs (SyncFinder::normalize_sync_quality (sync_quality));

          if (sync_quality > best_score.quality)
            {
              best_score.quality = sync_quality;
              best_score.speed = relative_speed * center;
            }
        }
    }
  std::lock_guard<std::mutex> lg (mutex);
  result_scores.push_back (best_score);
}

/*
 * The scores from speed search are usually a bit noisy, so the local maximum from the scores
 * vector is not necessarily the best choice.
 *
 * To get rid of the noise to some degree, this function smoothes the scores using a
 * cosine window and then finds the local maximum of this smooth function.
 */
static double
score_smooth_find_best (const vector<SpeedSync::Score>& in_scores, double step, double distance)
{
  auto scores = in_scores;
  sort (scores.begin(), scores.end(), [] (auto s1, auto s2) { return s1.speed < s2.speed; });

  double best_speed = 0;
  double best_quality = 0;

  for (double speed = scores.front().speed; speed < scores.back().speed; speed += 0.000001)
    {
      double quality_sum = 0;
      double quality_div = 0;

      for (auto s : scores)
        {
          double w = window_cos ((s.speed - speed) / (step * distance));

          quality_sum += s.quality * w;
          quality_div += w;
        }
      quality_sum /= quality_div;
      if (quality_sum > best_quality)
        {
          best_speed = speed;
          best_quality = quality_sum;
        }
    }

  return best_speed;
}

class SpeedSearch
{
  ThreadPool thread_pool;
  vector<std::unique_ptr<SpeedSync>> speed_sync;
  const WavData& in_data;
  double clip_location;

  double start_time;
  vector<double> times;

  SpeedSync *
  find_closest_speed_sync (double speed)
  {
    auto it = std::min_element (speed_sync.begin(), speed_sync.end(), [&](auto& x, auto& y)
      {
        return fabs (x->center_speed() - speed) < fabs (y->center_speed() - speed);
      });
    return (*it).get();
  }

  void
  timer_start()
  {
    start_time = get_time();
  }
  void
  timer_report()
  {
    auto t = get_time();
    times.push_back (t - start_time);
    start_time = t;
  }
public:
  SpeedSearch (const WavData& in_data, double clip_location) :
    in_data (in_data),
    clip_location (clip_location)
  {
  }
  vector<double>
  get_times()
  {
    return times;
  }
  static void
  debug_range (const SpeedScanParams& scan_params)
  {
    auto bound = [&] (float f) {
      return 100 * pow (scan_params.step, f * (scan_params.n_center_steps * (scan_params.n_steps * 2 + 1) + scan_params.n_steps));
    };
    printf ("range = [ %.2f .. %.2f ]\n", bound (-1), bound (1));
  }

  vector<SpeedSync::Score> run_search (const SpeedScanParams& scan_params, const vector<double>& speeds);
  vector<SpeedSync::Score> refine_search (const SpeedScanParams& scan_params, double speed);
};

vector<SpeedSync::Score>
SpeedSearch::run_search (const SpeedScanParams& scan_params, const vector<double>& speeds)
{
  /* speed is between 0.8 and 1.25, so we use a clip seconds factor of 1.3 to provide enough samples */
  WavData in_clip = get_speed_clip (clip_location, in_data, scan_params.seconds * 1.3);

  speed_sync.clear();

  for (auto speed : speeds)
    {
      for (int c = -scan_params.n_center_steps; c <= scan_params.n_center_steps; c++)
        {
          double c_speed = speed * pow (scan_params.step, c * (scan_params.n_steps * 2 + 1));

          speed_sync.push_back (std::make_unique<SpeedSync> (in_clip, c_speed));
        }
    }

  timer_start();

  for (auto& s : speed_sync)
    s->start_prepare_job (thread_pool, scan_params);
  thread_pool.wait_all();

  timer_report();

  for (auto& s : speed_sync)
    s->start_search_jobs (thread_pool, scan_params, s->center_speed());
  thread_pool.wait_all();

  timer_report();

  vector<SpeedSync::Score> scores;
  for (auto& s : speed_sync)
    {
      vector<SpeedSync::Score> step_scores = s->get_scores();
      scores.insert (scores.end(), step_scores.begin(), step_scores.end());
    }
  return scores;
}

vector<SpeedSync::Score>
SpeedSearch::refine_search (const SpeedScanParams& scan_params, double speed)
{
  SpeedSync *center_speed_sync = find_closest_speed_sync (speed);

  timer_start();

  center_speed_sync->start_search_jobs (thread_pool, scan_params, speed);
  thread_pool.wait_all();

  timer_report();

  return center_speed_sync->get_scores();
}


static void
select_n_best_scores (vector<SpeedSync::Score>& scores, size_t n)
{
  sort (scores.begin(), scores.end(), [](auto a, auto b) { return a.speed < b.speed; });

  auto get_quality = [&] (int pos) // handle corner cases
    {
      if (pos >= 0 && size_t (pos) < scores.size())
        return scores[pos].quality;
      else
        return 0.0;
    };

  vector<SpeedSync::Score> lmax_scores;
  for (int x = 0; size_t (x) < scores.size(); x++)
    {
      /* check for peaks
       *  - single peak : quality of the middle value is larger than the quality of the left and right neighbour
       *  - double peak : two values have equal quality, this must be larger than left and right neighbour
       */
      const double q1 = get_quality (x - 1);
      const double q2 = get_quality (x);
      const double q3 = get_quality (x + 1);

      if (q1 <= q2 && q2 >= q3)
        {
          lmax_scores.push_back (scores[x]);
          x++; // score with quality q3 cannot be a local maximum
        }
    }
  sort (lmax_scores.begin(), lmax_scores.end(), [](auto a, auto b) { return a.quality > b.quality; });

  if (lmax_scores.size() > n)
    lmax_scores.resize (n);
  scores = lmax_scores;
}

static vector<double>
get_clip_locations (const WavData& in_data, int n)
{
  Random rng (0, Random::Stream::speed_clip);

  /* to improve performance, we don't hash all samples but just a few */
  const vector<float>& samples = in_data.samples();
  vector<float> xsamples;
  for (size_t p = 0; p < samples.size(); p += rng() % 1000)
    xsamples.push_back (samples[p]);

  rng.seed (Random::seed_from_hash (xsamples), Random::Stream::speed_clip);

  /* return a set of n possible clip locations */
  vector<double> result;
  for (int c = 0; c < n; c++)
    result.push_back (rng.random_double());
  return result;
}

static double
get_best_clip_location (const WavData& in_data, double seconds, int candidates)
{
  double clip_location = 0;
  double best_energy = 0;

  /* try a few clip locations, use the one with highest signal energy */
  for (auto location : get_clip_locations (in_data, candidates))
    {
      WavData wd = get_speed_clip (location, in_data, seconds);

      double energy = 0;
      for (auto s : wd.samples())
        energy += s * s;
      if (energy > best_energy)
        {
          best_energy = energy;
          clip_location = location;
        }
    }
  return clip_location;
}

double
detect_speed (const WavData& in_data, bool print_results)
{
  /* typically even for high strength we need at least a few seconds of audio
   * in in_data for successful speed detection, but our algorithm won't work at
   * all for very short input files
   */
  double in_seconds = double (in_data.n_frames()) / in_data.sample_rate();
  if (in_seconds < 0.25)
    return 1;

  const SpeedScanParams scan1_normal /* first pass: find approximation: speed approximately 0.8..1.25 */
    {
      .seconds        = 25,
      .step           = 1.0007,
      .n_steps        = 5,
      .n_center_steps = 28,
    };
  const SpeedScanParams scan1_patient
    {
      .seconds        = 50,
      .step           = 1.00035,
      .n_steps        = 11,
      .n_center_steps = 28,
    };
  const SpeedScanParams scan1 = Params::detect_speed_patient ? scan1_patient : scan1_normal;

  const SpeedScanParams scan2_normal /* second pass: improve approximation */
    {
      .seconds        = 50,
      .step           = 1.00035,
      .n_steps        = 1,
    };
  const SpeedScanParams scan2_patient
    {
      .seconds        = 50,
      .step           = 1.000175,
      .n_steps        = 1,
    };
  const SpeedScanParams scan2 = Params::detect_speed_patient ? scan2_patient : scan2_normal;

  const SpeedScanParams scan3 /* third pass: fast refine (not always perfect) */
    {
      .seconds        = 50,
      .step           = 1.00005,
      .n_steps        = 40,
    };
  const double scan3_smooth_distance = 20;
  const double speed_sync_threshold = 0.4;
  const int    n_best = Params::detect_speed_patient ? 15 : 5;

  // SpeedSearch::debug_range (scan1);

  const int    clip_candidates = 5;
  const double clip_location = get_best_clip_location (in_data, scan1.seconds, clip_candidates);

  vector<SpeedSync::Score> scores;
  SpeedSearch speed_search (in_data, clip_location);

  /* initial search using grid */
  scores = speed_search.run_search (scan1, { 1.0 });

  /* improve N best matches */
  select_n_best_scores (scores, n_best);

  vector<double> speeds;
  for (auto score : scores)
    speeds.push_back (score.speed);

  scores = speed_search.run_search (scan2, speeds);

  /* improve or refine best match */
  select_n_best_scores (scores, 1);
  if (Params::detect_speed_patient)
    {
      // slower version: prepare magnitudes again, according to best speed
      scores = speed_search.run_search (scan3, { scores[0].speed });
    }
  else
    {
      // faster version: keep already computed magnitudes
      scores = speed_search.refine_search (scan3, scores[0].speed);
    }
  double best_speed = score_smooth_find_best (scores, 1 - scan3.step, scan3_smooth_distance);

  double best_quality = 0;
  for (auto score : scores)
    best_quality = max (best_quality, score.quality);

  if (print_results)
    {
      double delta = -1;
      if (Params::test_speed > 0)
        delta = 100 * fabs (best_speed - Params::test_speed) / Params::test_speed;
      printf ("detect_speed %f %f %.4f   ", best_speed, best_quality, delta);

      double total = 0.0;
      for (auto t : speed_search.get_times())
        {
          total += t;
          printf (" %.3f", t);
        }
      printf (" %.3f\n", total);
    }

  if (best_quality > speed_sync_threshold)
    return best_speed;
  else
    return 1;
}
