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
#include "wmcommon.hh"

using std::complex;
using std::vector;
using std::string;
using std::min;

void
SyncFinder::init_up_down (const WavData& wav_data, Mode mode)
{
  sync_bits.clear();

  // "long" blocks consist of two "normal" blocks, which means
  //   the sync bits pattern is repeated after the end of the first block
  const int first_block_end = mark_sync_frame_count() + mark_data_frame_count();
  const int block_count = mode == Mode::CLIP ? 2 : 1;
  size_t n_bands = Params::max_band - Params::min_band + 1;

  UpDownGen up_down_gen (Random::Stream::sync_up_down);
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
              frame_bit.frame = sync_frame_pos (f + bit * Params::sync_frames_per_bit) + block * first_block_end;
              for (int ch = 0; ch < wav_data.n_channels(); ch++)
                {
                  if (block == 0)
                    {
                      for (auto u : frame_up)
                        frame_bit.up.push_back (u - Params::min_band + n_bands * ch);
                      for (auto d : frame_down)
                        frame_bit.down.push_back (d - Params::min_band + n_bands * ch);
                    }
                  else
                    {
                      for (auto u : frame_up)
                        frame_bit.down.push_back (u - Params::min_band + n_bands * ch);
                      for (auto d : frame_down)
                        frame_bit.up.push_back (d - Params::min_band + n_bands * ch);
                    }
                }
              std::sort (frame_bit.up.begin(), frame_bit.up.end());
              std::sort (frame_bit.down.begin(), frame_bit.down.end());
              frame_bits.push_back (frame_bit);
            }
        }
      std::sort (frame_bits.begin(), frame_bits.end(), [] (FrameBit& f1, FrameBit& f2) { return f1.frame < f2.frame; });
      sync_bits.push_back (frame_bits);
    }
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
SyncFinder::sync_decode (const WavData& wav_data, const size_t start_frame,
                         const vector<float>& fft_out_db,
                         const vector<char>&  have_frames,
                         ConvBlockType *block_type)
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
              const int index = ((start_frame + frame_bit.frame) * wav_data.n_channels()) * n_bands;
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

vector<SyncFinder::Score>
SyncFinder::search_approx (const WavData& wav_data, Mode mode)
{
  vector<float> fft_db;
  vector<char>  have_frames;
  vector<Score> sync_scores;

  // compute multiple time-shifted fft vectors
  size_t n_bands = Params::max_band - Params::min_band + 1;
  int total_frame_count = mark_sync_frame_count() + mark_data_frame_count();
  if (mode == Mode::CLIP)
    total_frame_count *= 2;
  for (size_t sync_shift = 0; sync_shift < Params::frame_size; sync_shift += Params::sync_search_step)
    {
      sync_fft (wav_data, sync_shift, frame_count (wav_data) - 1, fft_db, have_frames, /* want all frames */ {});
      for (int start_frame = 0; start_frame < frame_count (wav_data); start_frame++)
        {
          const size_t sync_index = start_frame * Params::frame_size + sync_shift;
          if ((start_frame + total_frame_count) * wav_data.n_channels() * n_bands < fft_db.size())
            {
              ConvBlockType block_type;
              double quality = sync_decode (wav_data, start_frame, fft_db, have_frames, &block_type);
              // printf ("%zd %f\n", sync_index, quality);
              sync_scores.emplace_back (Score { sync_index, quality, block_type });
            }
        }
    }
  sort (sync_scores.begin(), sync_scores.end(), [] (const Score& a, const Score &b) { return a.index < b.index; });
  return sync_scores;
}

void
SyncFinder::sync_select_by_threshold (vector<Score>& sync_scores)
{
  /* for strength 8 and above:
   *   -> more false positive candidates are rejected, so we can use a lower threshold
   *
   * for strength 7 and below:
   *   -> we need a higher threshold, because otherwise watermark detection takes too long
   */
  const double strength = Params::water_delta * 1000;
  const double sync_threshold1 = strength > 7.5 ? 0.4 : 0.5;

  vector<Score> selected_scores;

  for (size_t i = 0; i < sync_scores.size(); i++)
    {
      if (sync_scores[i].quality > sync_threshold1)
        {
          double q_last = -1;
          double q_next = -1;
          if (i > 0)
            q_last = sync_scores[i - 1].quality;

          if (i + 1 < sync_scores.size())
            q_next = sync_scores[i + 1].quality;

          if (sync_scores[i].quality >= q_last && sync_scores[i].quality >= q_next)
            {
              selected_scores.emplace_back (sync_scores[i]);
              i++; // score with quality q_next cannot be a local maximum
            }
        }
    }
  sync_scores = selected_scores;
}

void
SyncFinder::sync_select_n_best (vector<Score>& sync_scores, size_t n)
{
  std::sort (sync_scores.begin(), sync_scores.end(), [](Score& s1, Score& s2) { return s1.quality > s2.quality; });
  if (sync_scores.size() > n)
    sync_scores.resize (n);
}

void
SyncFinder::search_refine (const WavData& wav_data, Mode mode, vector<Score>& sync_scores)
{
  vector<float> fft_db;
  vector<char>  have_frames;
  vector<Score> result_scores;

  int total_frame_count = mark_sync_frame_count() + mark_data_frame_count();
  const int first_block_end = total_frame_count;
  if (mode == Mode::CLIP)
    total_frame_count *= 2;

  vector<char> want_frames (total_frame_count);
  for (size_t f = 0; f < mark_sync_frame_count(); f++)
    {
      want_frames[sync_frame_pos (f)] = 1;
      if (mode == Mode::CLIP)
        want_frames[first_block_end + sync_frame_pos (f)] = 1;
    }

  for (const auto& score : sync_scores)
    {
      //printf ("%zd %s %f", score.index, find_closest_sync (score.index).c_str(), score.quality);

      // refine match
      double best_quality       = score.quality;
      size_t best_index         = score.index;
      ConvBlockType best_block_type = score.block_type; /* doesn't really change during refinement */

      int start = std::max (int (score.index) - Params::sync_search_step, 0);
      int end   = score.index + Params::sync_search_step;
      for (int fine_index = start; fine_index <= end; fine_index += Params::sync_search_fine)
        {
          sync_fft (wav_data, fine_index, total_frame_count, fft_db, have_frames, want_frames);
          if (fft_db.size())
            {
              ConvBlockType block_type;
              double        q = sync_decode (wav_data, 0, fft_db, have_frames, &block_type);

              if (q > best_quality)
                {
                  best_quality = q;
                  best_index   = fine_index;
                }
            }
        }
      //printf (" => refined: %zd %s %f\n", best_index, find_closest_sync (best_index).c_str(), best_quality);
      if (best_quality > Params::sync_threshold2)
        result_scores.push_back (Score { best_index, best_quality, best_block_type });
    }
  sync_scores = result_scores;
}

vector<SyncFinder::Score>
SyncFinder::fake_sync (const WavData& wav_data, Mode mode)
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

  return result_scores;
}

vector<SyncFinder::Score>
SyncFinder::search (const WavData& wav_data, Mode mode)
{
  if (Params::test_no_sync)
    return fake_sync (wav_data, mode);

  init_up_down (wav_data, mode);

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
  vector<Score> sync_scores = search_approx (wav_data, mode);

  sync_select_by_threshold (sync_scores);
  if (mode == Mode::CLIP)
    sync_select_n_best (sync_scores, 5);

  search_refine (wav_data, mode, sync_scores);

  return sync_scores;
}

vector<vector<SyncFinder::FrameBit>>
SyncFinder::get_sync_bits (const WavData& wav_data, Mode mode)
{
  init_up_down (wav_data, mode);
  return sync_bits;
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

  fft_out_db.resize (wav_data.n_channels() * n_bands * frame_count);
  have_frames.resize (frame_count);

  for (size_t f = 0; f < frame_count; f++)
    {
      const size_t f_first = (index + f * Params::frame_size) * wav_data.n_channels();
      const size_t f_last  = (index + (f + 1) * Params::frame_size) * wav_data.n_channels();

      if ((want_frames.size() && !want_frames[f])   // frame not wanted?
      ||  (f_last < wav_data_first)                 // frame in silence before input?
      ||  (f_first > wav_data_last))                // frame in silence after input?
        {
          out_pos += n_bands * wav_data.n_channels();
        }
      else
        {
          constexpr double min_db = -96;

          vector<vector<complex<float>>> frame_result = fft_analyzer.run_fft (samples, index + f * Params::frame_size);

          /* computing db-magnitude is expensive, so we better do it here */
          for (int ch = 0; ch < wav_data.n_channels(); ch++)
            for (int i = Params::min_band; i <= Params::max_band; i++)
              fft_out_db[out_pos++] = db_from_complex (frame_result[ch][i], min_db);

          have_frames[f] = 1;
        }
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
