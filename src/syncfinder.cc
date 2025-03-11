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

#include <vector>
#include <algorithm>

#include "syncfinder.hh"
#include "threadpool.hh"
#include "wmcommon.hh"

using std::complex;
using std::vector;
using std::string;
using std::min;

vector<vector<SyncFinder::FrameBit>>
SyncFinder::get_sync_bits (const Key& key, Mode mode)
{
  vector<vector<SyncFinder::FrameBit>> sync_bits;

  // "long" blocks consist of two "normal" blocks, which means
  //   the sync bits pattern is repeated after the end of the first block
  const int first_block_end = mark_sync_frame_count() + mark_data_frame_count();
  const int block_count = mode == Mode::CLIP ? 2 : 1;

  UpDownGen up_down_gen (key, Random::Stream::sync_up_down);
  BitPosGen bit_pos_gen (key);
  for (int bit = 0; bit < Params::sync_bits; bit++)
    {
      vector<FrameBit> frame_bits;
      for (int f = 0; f < Params::sync_frames_per_bit; f++)
        {
          UpDownArray frame_up, frame_down;
          up_down_gen.get (f + bit * Params::sync_frames_per_bit, frame_up, frame_down);

          for (int block = 0; block < block_count; block++)
            {
              FrameBit frame_bit;
              frame_bit.frame = bit_pos_gen.sync_frame (f + bit * Params::sync_frames_per_bit) + block * first_block_end;
              if (block == 0)
                {
                  for (auto u : frame_up)
                    frame_bit.up.push_back (u - Params::min_band);
                  for (auto d : frame_down)
                    frame_bit.down.push_back (d - Params::min_band);
                }
              else
                {
                  for (auto u : frame_up)
                    frame_bit.down.push_back (u - Params::min_band);
                  for (auto d : frame_down)
                    frame_bit.up.push_back (d - Params::min_band);
                }
              std::sort (frame_bit.up.begin(), frame_bit.up.end());
              std::sort (frame_bit.down.begin(), frame_bit.down.end());
              frame_bits.push_back (frame_bit);
            }
        }
      std::sort (frame_bits.begin(), frame_bits.end(), [] (FrameBit& f1, FrameBit& f2) { return f1.frame < f2.frame; });
      sync_bits.push_back (frame_bits);
    }
  return sync_bits;
}

/* safe to call from any thread */
double
SyncFinder::normalize_sync_quality (double raw_quality)
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

/* safe to call from any thread */
double
SyncFinder::bit_quality (float umag, float dmag, int bit)
{
  const int expect_data_bit = bit & 1; /* expect 010101 */

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
  return expect_data_bit ? raw_bit : -raw_bit;
}

double
SyncFinder::sync_decode (const vector<vector<FrameBit>>& sync_bits,
                         const size_t start_frame,
                         const vector<float>& fft_out_db,
                         const vector<char>&  have_frames)
{
  double sync_quality = 0;

  size_t n_bands = Params::max_band - Params::min_band + 1;
  int bit_count = 0;
  for (size_t bit = 0; bit < sync_bits.size(); bit++)
    {
      const vector<FrameBit>& frame_bits = sync_bits[bit];
      float umag = 0, dmag = 0;

      int frame_bit_count = 0;
      for (const auto& frame_bit : frame_bits)
        {
          if (have_frames[start_frame + frame_bit.frame])
            {
              const int index = (start_frame + frame_bit.frame) * n_bands;
              for (size_t i = 0; i < frame_bit.up.size(); i++)
                {
                  umag += fft_out_db[index + frame_bit.up[i]];
                  dmag += fft_out_db[index + frame_bit.down[i]];
                }
              frame_bit_count++;
            }
        }
      sync_quality += bit_quality (umag, dmag, bit) * frame_bit_count;
      bit_count += frame_bit_count;
    }
  if (bit_count)
    sync_quality /= bit_count;
  sync_quality = normalize_sync_quality (sync_quality);

  return sync_quality;
}

void
SyncFinder::scan_silence (const WavData& wav_data)
{
  const vector<float>& samples = wav_data.samples();

  // find first non-zero sample
  wav_data_first = 0;
  while (wav_data_first < samples.size() && samples[wav_data_first] == 0)
    wav_data_first++;

  // search wav_data_last to get [wav_data_first, wav_data_last) range
  wav_data_last = samples.size();
  while (wav_data_last > wav_data_first && samples[wav_data_last - 1] == 0)
    wav_data_last--;
}

void
SyncFinder::search_approx (vector<SearchKeyResult>& key_results, const vector<vector<vector<FrameBit>>>& sync_bits, const WavData& wav_data, Mode mode)
{
  ThreadPool    thread_pool;
  vector<float> fft_db;
  vector<char>  have_frames;

  std::mutex result_mutex;
  // compute multiple time-shifted fft vectors
  size_t n_bands = Params::max_band - Params::min_band + 1;
  int total_frame_count = mark_sync_frame_count() + mark_data_frame_count();
  if (mode == Mode::CLIP)
    total_frame_count *= 2;
  for (size_t sync_shift = 0; sync_shift < Params::frame_size; sync_shift += Params::sync_search_step)
    {
      sync_fft_parallel (thread_pool, wav_data, sync_shift, fft_db, have_frames);

      vector<int> start_frames;
      for (int start_frame = 0; start_frame < frame_count (wav_data); start_frame++)
        {
          if ((start_frame + total_frame_count) * n_bands < fft_db.size())
            start_frames.push_back (start_frame);
        }

      for (size_t k = 0; k < key_results.size(); k++)
        {
          for (auto split_start_frames : split_vector (start_frames, 256))
            {
              thread_pool.add_job ([this, k, sync_shift, split_start_frames,
                                    &sync_bits, &wav_data, &fft_db, &have_frames, &key_results, &result_mutex]()
                {
                  for (auto start_frame : split_start_frames)
                    {
                      double quality = sync_decode (sync_bits[k], start_frame, fft_db, have_frames);
                      // printf ("%zd %f\n", sync_index, quality);
                      const size_t sync_index = start_frame * Params::frame_size + sync_shift;
                      {
                        std::lock_guard<std::mutex> lg (result_mutex);
                        SearchScore search_score;
                        search_score.index       = sync_index;
                        search_score.raw_quality = quality;
                        search_score.local_mean  = 0; // fill this after all search scores are ready
                        key_results[k].scores.push_back (search_score);
                      }
                    }
                });
            }
        }
      thread_pool.wait_all();
    }
  for (auto& key_result : key_results)
    {
      sort (key_result.scores.begin(), key_result.scores.end(), [] (const SearchScore& a, const SearchScore &b) { return a.index < b.index; });

      /* compute local mean for all scores */
      for (int i = 0; i < int (key_result.scores.size()); i++)
        {
          double avg = 0;
          int n = 0;
          for (int j = -20; j <= 20; j++)
            {
              if (std::abs (j) >= 4)
                {
                  int idx = i + j;
                  if (idx >= 0 && idx < int (key_result.scores.size()))
                    {
                      avg += key_result.scores[idx].raw_quality;
                      n++;
                    }
                }
            }
          if (n > 0)
            avg /= n;
          if (getenv ("NOAVG"))
            avg = 0;
          key_result.scores[i].local_mean = avg;
          //printf ("%f %f #Q\n", key_result.scores[i].raw_quality, key_result.scores[i].local_mean);
        }
    }
}

void
SyncFinder::sync_select_by_threshold (vector<SearchScore>& sync_scores)
{
  const double sync_threshold1 = Params::sync_threshold2 * 0.75;

  vector<SearchScore> selected_scores;

  for (size_t i = 0; i < sync_scores.size(); i++)
    {
      double q = sync_scores[i].abs_quality();
      if (q > sync_threshold1)
        {
          double q_last = 0;
          double q_next = 0;
          if (i > 0)
            q_last = sync_scores[i - 1].abs_quality();

          if (i + 1 < sync_scores.size())
            q_next = sync_scores[i + 1].abs_quality();

          if (q >= q_last && q >= q_next)
            {
              selected_scores.emplace_back (sync_scores[i]);
              i++; // score with quality q_next cannot be a local maximum
            }
        }
    }
  sync_scores = selected_scores;
}

void
SyncFinder::sync_select_n_best (vector<SearchScore>& sync_scores, size_t n)
{
  std::sort (sync_scores.begin(), sync_scores.end(), [](SearchScore& s1, SearchScore& s2) { return s1.abs_quality() > s2.abs_quality(); });
  if (sync_scores.size() > n)
    sync_scores.resize (n);
}

void
SyncFinder::search_refine (const WavData& wav_data, Mode mode, SearchKeyResult& key_result, const vector<vector<FrameBit>>& sync_bits)
{
  ThreadPool          thread_pool;
  std::mutex          result_mutex;
  vector<SearchScore> result_scores;
  BitPosGen           bit_pos_gen (key_result.key);

  int total_frame_count = mark_sync_frame_count() + mark_data_frame_count();
  const int first_block_end = total_frame_count;
  if (mode == Mode::CLIP)
    total_frame_count *= 2;

  vector<char> want_frames (total_frame_count);
  for (size_t f = 0; f < mark_sync_frame_count(); f++)
    {
      want_frames[bit_pos_gen.sync_frame (f)] = 1;
      if (mode == Mode::CLIP)
        want_frames[first_block_end + bit_pos_gen.sync_frame (f)] = 1;
    }

  for (const auto& score : key_result.scores)
    {
      thread_pool.add_job ([this, score, total_frame_count,
                            &wav_data, &want_frames, &sync_bits, &result_scores, &result_mutex] ()
        {
          vector<float> fft_db;
          vector<char>  have_frames;
          //printf ("%zd %s %f", score.index, find_closest_sync (score.index).c_str(), score.quality);

          // refine match
          double best_quality       = 0;
          size_t best_index         = score.index;

          int start = std::max (int (score.index) - Params::sync_search_step, 0);
          int end   = score.index + Params::sync_search_step;
          for (int fine_index = start; fine_index <= end; fine_index += Params::sync_search_fine)
            {
              sync_fft (wav_data, fine_index, total_frame_count, fft_db, have_frames, want_frames);
              if (fft_db.size())
                {
                  double q = sync_decode (sync_bits, 0, fft_db, have_frames);

                  if (fabs (q) > fabs (best_quality))
                    {
                      best_quality = q;
                      best_index   = fine_index;
                    }
                }
            }
          //printf (" => refined: %zd %s %f\n", best_index, find_closest_sync (best_index).c_str(), best_quality);
          if (fabs (best_quality - score.local_mean) > Params::sync_threshold2)
            {
              std::lock_guard<std::mutex> lg (result_mutex);
              SearchScore refined_score;
              refined_score.index = best_index;
              refined_score.raw_quality = best_quality;
              refined_score.local_mean = score.local_mean;

              result_scores.push_back (refined_score);
            }
        });
    }
  thread_pool.wait_all();
  sort (result_scores.begin(), result_scores.end(), [] (const SearchScore& a, const SearchScore &b) { return a.index < b.index; });
  key_result.scores = result_scores;
}

vector<SyncFinder::KeyResult>
SyncFinder::fake_sync (const vector<Key>& key_list, const WavData& wav_data, Mode mode)
{
  vector<Score> result_scores;

  if (mode == Mode::BLOCK)
    {
      const size_t expect0 = Params::frames_pad_start * Params::frame_size;
      const size_t expect_step = (mark_sync_frame_count() + mark_data_frame_count()) * Params::frame_size;
      const size_t expect_end = frame_count (wav_data) * Params::frame_size;

      int ab = 0;
      for (size_t expect_index = expect0; expect_index + expect_step < expect_end; expect_index += expect_step)
        result_scores.push_back (Score { expect_index, 1.0, (ab++ & 1) ? ConvBlockType::b : ConvBlockType::a });
    }

  vector<KeyResult> key_results;
  for (auto key : key_list)
    {
      KeyResult key_result;
      key_result.key = key;
      key_result.sync_scores = result_scores;
      key_results.push_back (key_result);
    }
  return key_results;
}

vector<SyncFinder::KeyResult>
SyncFinder::search (const vector<Key>& key_list, const WavData& wav_data, Mode mode)
{
  if (Params::test_no_sync)
    return fake_sync (key_list, wav_data, mode);

  if (mode == Mode::CLIP)
    {
      /* in clip mode we optimize handling large areas of padding which is silent */
      scan_silence (wav_data);
    }
  else
    {
      /* in block mode we don't do anything special for silence at beginning/end */
      wav_data_first = 0;
      wav_data_last  = wav_data.samples().size();
    }

  vector<SearchKeyResult>           search_key_results;
  vector<vector<vector<FrameBit>>>  sync_bits;

  for (const auto& key : key_list)
    {
      SearchKeyResult search_key_result;
      search_key_result.key = key;
      search_key_results.push_back (search_key_result);
      sync_bits.push_back (get_sync_bits (key, mode));
    }

  search_approx (search_key_results, sync_bits, wav_data, mode);
  vector<SyncFinder::KeyResult> key_results;
  for (size_t k = 0; k < search_key_results.size(); k++)
    {
      /* find local maxima, select by threshold */
      sync_select_by_threshold (search_key_results[k].scores);
      if (mode == Mode::CLIP)
        sync_select_n_best (search_key_results[k].scores, 5);

      search_refine (wav_data, mode, search_key_results[k], sync_bits[k]);

      KeyResult key_result;
      key_result.key = search_key_results[k].key;
      for (auto search_score : search_key_results[k].scores)
        {
          double q = search_score.raw_quality - search_score.local_mean;

          Score score;
          score.index = search_score.index;
          score.quality = fabs (q);
          score.block_type = q > 0 ? ConvBlockType::a : ConvBlockType::b;
          key_result.sync_scores.push_back (score);
        }
      key_results.push_back (key_result);
    }

  return key_results;
}

void
SyncFinder::sync_fft (const WavData& wav_data, size_t index, size_t frame_count, vector<float>& fft_out_db, vector<char>& have_frames, const vector<char>& want_frames)
{
  fft_out_db.clear();
  have_frames.clear();

  /* read past end? -> fail */
  if (wav_data.n_values() < (index + frame_count * Params::frame_size) * wav_data.n_channels())
    return;

  FFTAnalyzer fft_analyzer (wav_data.n_channels());
  const vector<float>& samples = wav_data.samples();
  const size_t n_bands = Params::max_band - Params::min_band + 1;
  int out_pos = 0;

  fft_out_db.resize (n_bands * frame_count);
  have_frames.resize (frame_count);

  for (size_t f = 0; f < frame_count; f++)
    {
      const size_t f_first = (index + f * Params::frame_size) * wav_data.n_channels();
      const size_t f_last  = (index + (f + 1) * Params::frame_size) * wav_data.n_channels();

      if ((want_frames.size() && !want_frames[f])   // frame not wanted?
      ||  (f_last < wav_data_first)                 // frame in silence before input?
      ||  (f_first > wav_data_last))                // frame in silence after input?
        {
          out_pos += n_bands;
        }
      else
        {
          constexpr double min_db = -96;

          vector<vector<complex<float>>> frame_result = fft_analyzer.run_fft (samples, index + f * Params::frame_size);

          /* computing db-magnitude is expensive, so we better do it here */
          for (int ch = 0; ch < wav_data.n_channels(); ch++)
            for (int i = Params::min_band; i <= Params::max_band; i++)
              fft_out_db[out_pos + i - Params::min_band] += db_from_complex (frame_result[ch][i], min_db);

          out_pos += n_bands;

          have_frames[f] = 1;
        }
    }
}

void
SyncFinder::sync_fft_parallel (ThreadPool& thread_pool,
                               const WavData& wav_data,
                               size_t index,
                               std::vector<float>& fft_out_db,
                               std::vector<char>& have_frames)
{
  std::mutex result_mutex;

  fft_out_db.clear();
  have_frames.clear();

  struct PartialFFTResult
  {
    int           start_frame = 0;
    vector<char>  have_frames;
    vector<float> fft_db;
  };
  vector<PartialFFTResult> partial_fft_results;
  const int frames_per_job = 256;
  for (int start_frame = 0; start_frame < frame_count (wav_data); start_frame += frames_per_job)
    {
      thread_pool.add_job ([this, start_frame, index, frames_per_job,
                            &wav_data, &partial_fft_results, &result_mutex]
        {
          const int remaining_frames = frame_count (wav_data) - 1 - start_frame;
          const int frames = std::min (remaining_frames, frames_per_job);
          if (frames > 0)
            {
              PartialFFTResult result;
              result.start_frame = start_frame;
              sync_fft (wav_data, index + start_frame * Params::frame_size, frames, result.fft_db, result.have_frames, /* want all frames */ {});
              if (!result.fft_db.size())
                warning ("SyncFinder: sync_fft_parallel expected %d fft frames, but result was empty\n", frames);
              {
                std::lock_guard<std::mutex> lg (result_mutex);
                partial_fft_results.push_back (result);
              }
            }
        });
    }
  thread_pool.wait_all();
  std::sort (partial_fft_results.begin(), partial_fft_results.end(),
             [] (const auto& r1, const auto& r2) { return r1.start_frame < r2.start_frame; });

  for (const auto& partial_fft_result : partial_fft_results)
    {
      fft_out_db.insert (fft_out_db.end(), partial_fft_result.fft_db.begin(), partial_fft_result.fft_db.end());
      have_frames.insert (have_frames.end(), partial_fft_result.have_frames.begin(), partial_fft_result.have_frames.end());
    }
}

string
SyncFinder::find_closest_sync (size_t index)
{
  int wm_length = (mark_data_frame_count() + mark_sync_frame_count()) * Params::frame_size;
  int wm_offset = Params::frames_pad_start * Params::frame_size;
  int best_error = wm_length * 2;
  int best = 0;

  for (int i = 0; i < 100; i++)
    {
      int error = abs (int (index) - (wm_offset + i * wm_length));
      if (error < best_error)
        {
          best = i;
          best_error = error;
        }
    }
  return string_printf ("n:%d offset:%d", best, int (index) - (wm_offset + best * wm_length));
}

vector<vector<int>>
SyncFinder::split_vector (vector<int>& in_vector, size_t max_size)
{
  /* split input vector into smaller vectors of at most max_size elements */
  vector<vector<int>> result;
  size_t size = in_vector.size();

  for (size_t i = 0; i < size; i += max_size)
    result.emplace_back (in_vector.begin() + i, in_vector.begin() + i + min (size - i, max_size));
  return result;
}
