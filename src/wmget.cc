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

#include <string>
#include <algorithm>

#include <zita-resampler/resampler.h>
#include <zita-resampler/vresampler.h>

#include "wavdata.hh"
#include "wmcommon.hh"
#include "convcode.hh"
#include "shortcode.hh"
#include "syncfinder.hh"
#include "fft.hh"

using std::string;
using std::vector;
using std::min;
using std::max;
using std::complex;

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
  error ("audiowmark: resampling from rate %d to rate %d not supported.\n", wav_data.sample_rate(), rate);
  exit (1);
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

static vector<float>
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

static vector<float>
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

static vector<float>
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

class ResultSet
{
public:
  enum class Type { BLOCK, CLIP, ALL };
  struct Pattern
  {
    vector<int>       bit_vec;
    float             decode_error = 0;
    SyncFinder::Score sync_score;
    Type              type;
  };
private:
  vector<Pattern> patterns;

public:
  void
  add_pattern (SyncFinder::Score sync_score, const vector<int>& bit_vec, float decode_error, Type pattern_type)
  {
    Pattern p;
    p.sync_score = sync_score;
    p.bit_vec = bit_vec;
    p.decode_error = decode_error;
    p.type = pattern_type;

    patterns.push_back (p);
  }
  void
  print()
  {
    std::stable_sort (patterns.begin(), patterns.end(), [](const Pattern& p1, const Pattern& p2) {
      const int all1 = p1.type == Type::ALL;
      const int all2 = p2.type == Type::ALL;
      if (all1 != all2)
        return all1 < all2;
      else
        return p1.sync_score.index < p2.sync_score.index;
    });
    for (const auto& pattern : patterns)
      {
        if (pattern.type == Type::ALL) /* this is the combined pattern "all" */
          {
            printf ("pattern   all %s %.3f %.3f\n", bit_vec_to_str (pattern.bit_vec).c_str(),
                                                    pattern.sync_score.quality, pattern.decode_error);
          }
        else
          {
            string block_str;

            switch (pattern.sync_score.block_type)
              {
                case ConvBlockType::a:  block_str = "A";
                                        break;
                case ConvBlockType::b:  block_str = "B";
                                        break;
                case ConvBlockType::ab: block_str = "AB";
                                        break;
              }
            if (pattern.type == Type::CLIP)
              block_str = "CLIP-" + block_str;

            const int seconds = pattern.sync_score.index / Params::mark_sample_rate;
            printf ("pattern %2d:%02d %s %.3f %.3f %s\n", seconds / 60, seconds % 60, bit_vec_to_str (pattern.bit_vec).c_str(),
                      pattern.sync_score.quality, pattern.decode_error, block_str.c_str());
          }
      }
  }
  int
  print_match_count (const string& orig_pattern)
  {
    int match_count = 0;

    vector<int> orig_vec = bit_str_to_vec (orig_pattern);
    for (auto p : patterns)
      {
        bool        match = true;

        for (size_t i = 0; i < p.bit_vec.size(); i++)
          match = match && (p.bit_vec[i] == orig_vec[i % orig_vec.size()]);

        if (match)
          match_count++;
      }
    printf ("match_count %d %zd\n", match_count, patterns.size());
    return match_count;
  }
  double
  best_quality() const
  {
    double q = -1;
    for (const auto& pattern : patterns)
      if (pattern.sync_score.quality > q)
        q = pattern.sync_score.quality;
    return q;
  }
};

/*
 * The block decoder is responsible for finding whole data blocks inside the
 * input file and decoding them. This only works for files that are large
 * enough to contain one data block. Incomplete blocks are ignored.
 *
 * INPUT:   AA|BBBBB|AAAAA|BBB
 * MATCH       BBBBB
 * MATCH             AAAAA
 *
 * The basic algorithm is this:
 *
 *  - use sync finder to find start index for the blocks
 *  - decode the blocks
 *  - try to combine A + B blocks for better error correction (AB)
 *  - try to combine all available blocks for better error correction (all pattern)
 */
class BlockDecoder
{
  int debug_sync_frame_count = 0;
  vector<SyncFinder::Score> sync_scores; // stored here for sync debugging
public:
  void
  run (const WavData& wav_data, ResultSet& result_set)
  {
    int total_count = 0;

    SyncFinder sync_finder;
    sync_scores = sync_finder.search (wav_data, SyncFinder::Mode::BLOCK);

    vector<float> raw_bit_vec_all (code_size (ConvBlockType::ab, Params::payload_size));
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
            assert (raw_bit_vec.size() == code_size (ConvBlockType::a, Params::payload_size));

            raw_bit_vec = randomize_bit_order (raw_bit_vec, /* encode */ false);

            /* ---- deal with this pattern ---- */
            float decode_error = 0;
            vector<int> bit_vec = code_decode_soft (sync_score.block_type, normalize_soft_bits (raw_bit_vec), &decode_error);

            if (!bit_vec.empty())
              result_set.add_pattern (sync_score, bit_vec, decode_error, ResultSet::Type::BLOCK);
            total_count += 1;

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
                vector<int> bit_vec = code_decode_soft (ConvBlockType::ab, normalize_soft_bits (ab_bits), &decode_error);
                if (!bit_vec.empty())
                  {
                    score_ab.index = sync_score.index;
                    score_ab.quality = (ab_quality[0] + ab_quality[1]) / 2;
                    result_set.add_pattern (score_ab, bit_vec, decode_error, ResultSet::Type::BLOCK);
                  }
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
        vector<int> bit_vec = code_decode_soft (ConvBlockType::ab, soft_bit_vec, &decode_error);

        if (!bit_vec.empty())
          result_set.add_pattern (score_all, bit_vec, decode_error, ResultSet::Type::ALL);
      }

    debug_sync_frame_count = frame_count (wav_data);
  }
  void
  print_debug_sync()
  {
      /* search sync markers at typical positions */
      const int expect0 = Params::frames_pad_start * Params::frame_size;
      const int expect_step = (mark_sync_frame_count() + mark_data_frame_count()) * Params::frame_size;
      const int expect_end = debug_sync_frame_count * Params::frame_size;

      int sync_match = 0;
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
};

/*
 * The clip decoder is responsible for decoding short clips. It is designed to
 * handle input sizes that are smaller than one data block. One case is that
 * the clip contains a partial A block (so the data could start after the start
 * of the A block and end before the end of the A block).
 *
 * ORIG:   |AAAAA|BBBBB|AAAAA|BBBBB|
 * CLIP:    |AAA|
 *
 * A clip could also contain the end of one block and the start of the next block,
 * like this:
 *
 * ORIG:   |AAAAA|BBBBB|AAAAA|BBBBB|
 * CLIP:                   |A|BB|
 *
 * The basic algorithm is this:
 *
 *  - zeropad   |AAA|  => 00000|AAA|00000
 *  - use sync finder to find start index of one long block in the zeropadded data
 *  - decode the bits
 *
 * For files larger than one data block, we decode twice, at the beginning and end
 *
 * INPUT   AAA|BBBBB|A
 * CLIP #1 AAA|BB
 * CLIP #2      BBBB|A
 */
class ClipDecoder
{
  const int frames_per_block = 0;

  vector<float>
  mix_or_linear_decode (vector<vector<complex<float>>>& fft_out, int n_channels)
  {
    if (Params::mix)
      return mix_decode (fft_out, n_channels);
    else
      return linear_decode (fft_out, n_channels);
  }
  void
  run_padded (const WavData& wav_data, ResultSet& result_set, double time_offset_sec)
  {
    SyncFinder                sync_finder;
    vector<SyncFinder::Score> sync_scores = sync_finder.search (wav_data, SyncFinder::Mode::CLIP);
    FFTAnalyzer               fft_analyzer (wav_data.n_channels());

    for (auto sync_score : sync_scores)
      {
        const size_t count = mark_sync_frame_count() + mark_data_frame_count();
        const size_t index = sync_score.index;
        auto fft_range_out1 = fft_analyzer.fft_range (wav_data.samples(), index, count);
        auto fft_range_out2 = fft_analyzer.fft_range (wav_data.samples(), index + count * Params::frame_size, count);
        if (fft_range_out1.size() && fft_range_out2.size())
          {
            const auto raw_bit_vec1 = randomize_bit_order (mix_or_linear_decode (fft_range_out1, wav_data.n_channels()), /* encode */ false);
            const auto raw_bit_vec2 = randomize_bit_order (mix_or_linear_decode (fft_range_out2, wav_data.n_channels()), /* encode */ false);
            const size_t bits_per_block = raw_bit_vec1.size();
            vector<float> raw_bit_vec;
            for (size_t i = 0; i < bits_per_block; i++)
              {
                if (sync_score.block_type == ConvBlockType::a)
                  {
                    raw_bit_vec.push_back (raw_bit_vec1[i]);
                    raw_bit_vec.push_back (raw_bit_vec2[i]);
                  }
                else
                  {
                    raw_bit_vec.push_back (raw_bit_vec2[i]);
                    raw_bit_vec.push_back (raw_bit_vec1[i]);
                  }
              }

            float decode_error = 0;
            vector<int> bit_vec = code_decode_soft (ConvBlockType::ab, normalize_soft_bits (raw_bit_vec), &decode_error);
            if (!bit_vec.empty())
              {
                SyncFinder::Score sync_score_nopad = sync_score;
                sync_score_nopad.index = time_offset_sec * wav_data.sample_rate();
                result_set.add_pattern (sync_score_nopad, bit_vec, decode_error, ResultSet::Type::CLIP);
              }
          }
      }
  }
  enum class Pos { START, END };
  void
  run_block (const WavData& wav_data, ResultSet& result_set, Pos pos)
  {
    const size_t n = (frames_per_block + 5) * Params::frame_size * wav_data.n_channels();

    // range of samples used by clip: [first_sample, last_sample)
    size_t first_sample;
    size_t last_sample;
    size_t pad_samples_start = n;
    size_t pad_samples_end   = n;

    if (pos == Pos::START)
      {
        first_sample = 0;
        last_sample  = min (n, wav_data.n_values());

        // increase padding at start for small blocks
        //   -> (available samples + padding) must always be one L-block
        if (last_sample < n)
          pad_samples_start += n - last_sample;
      }
    else // (pos == Pos::END)
      {
        if (wav_data.n_values() <= n)
          return;

        first_sample = wav_data.n_values() - n;
        last_sample  = wav_data.n_values();
      }
    const double time_offset = double (first_sample) / wav_data.sample_rate() / wav_data.n_channels();
    vector<float> ext_samples (wav_data.samples().begin() + first_sample, wav_data.samples().begin() + last_sample);

    if (0)
      {
        printf ("%d: %f..%f\n", int (pos), time_offset, time_offset + double (ext_samples.size()) / wav_data.sample_rate() / wav_data.n_channels());
        printf ("%f< >%f\n",
          double (pad_samples_start) / wav_data.sample_rate() / wav_data.n_channels(),
          double (pad_samples_end) / wav_data.sample_rate() / wav_data.n_channels());
      }
    ext_samples.insert (ext_samples.begin(), pad_samples_start, 0);
    ext_samples.insert (ext_samples.end(),   pad_samples_end, 0);

    WavData l_wav_data (ext_samples, wav_data.n_channels(), wav_data.sample_rate(), wav_data.bit_depth());
    run_padded (l_wav_data, result_set, time_offset);
   }
public:
  ClipDecoder() :
    frames_per_block (mark_sync_frame_count() + mark_data_frame_count())
  {
  }
  void
  run (const WavData& wav_data, ResultSet& result_set)
  {
    const int wav_frames = wav_data.n_values() / (Params::frame_size * wav_data.n_channels());
    if (wav_frames < frames_per_block * 3.1) /* clip decoder is only used for small wavs */
      {
        run_block (wav_data, result_set, Pos::START);
        run_block (wav_data, result_set, Pos::END);
      }
  }
};

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

static double
detect_speed (const WavData& wav_data, double center, double step, int n_steps, int seconds, double *quality_p)
{
  WavData wd_truncated = truncate (wav_data, seconds * 1.5);
  double best_speed = 1.0;
  double best_quality = 0;

  int best_hint_step = 0;
  if (Params::detect_speed_hint > 0)
    {
      double best_dist = 1000;
      for (int p = -n_steps; p <= n_steps; p++)
        {
          double dist = fabs (center * pow (step, p) - Params::detect_speed_hint);
          if (dist < best_dist)
            {
              best_hint_step = p;
              best_dist      = dist;
            }
        }
    }
  printf ("## range [%f..%f], n_steps=%d\n", center * pow (step, -n_steps), center * pow (step, n_steps), n_steps);
  for (int p = -n_steps; p <= n_steps; p++)
    {
      if (Params::detect_speed_hint > 0)
        if (abs (p - best_hint_step) > 2)
          continue;

      double speed = center * pow (step, p);

      WavData wd_resampled = resample_ratio (wd_truncated, speed, Params::mark_sample_rate);
      wd_resampled = truncate (wd_resampled, seconds);

      ResultSet result_set;
      ClipDecoder clip_decoder;
      clip_decoder.run (wd_resampled, result_set);
      printf ("%f %f         ", speed, result_set.best_quality());
      if (result_set.best_quality() > 0)
        {
          printf ("\n");
          if (result_set.best_quality() > best_quality)
            {
              best_quality = result_set.best_quality();
              best_speed = speed;
            }
        }
      else
        printf ("\r");
      fflush (stdout);
    }
  if (quality_p)
    *quality_p = best_quality;
  return best_speed;
}

class SpeedSync
{
public:
  struct Score
  {
    double speed = 0;;
    double quality = 0;
  };
private:
  struct Mags {
    float umag = 0;
    float dmag = 0;
  };
  vector<vector<SyncFinder::FrameBit>> sync_bits;
  vector<vector<Mags>> fft_sync_bits;
  void  prepare_mags (const WavData& in_data, double center, double seconds);
  Score compare (double relative_speed, double center);
public:
  vector<Score>
  search (const WavData& in_data, double center, double step, int n_steps, double seconds)
  {
    prepare_mags (in_data, center, seconds);

    vector<Score> scores;
    for (int p = -n_steps; p <= n_steps; p++)
      {
        const double relative_speed = pow (step, p);

        scores.push_back (compare (relative_speed, center));
      }
    return scores;
  }
};

static double
speed_scan (const WavData& in_data)
{
  vector<SpeedSync::Score> scores;

  /* n_center_steps / n_steps settings: speed approximately 0.8..1.25 */
  const int n_center_steps = 28;
  const int n_steps = 5;
  const double step = 1.0007;
  for (int c = -n_center_steps; c <= n_center_steps; c++)
    {
      double c_speed = pow (step, c * (n_steps * 2 + 1));

      SpeedSync speed_sync;
      vector<SpeedSync::Score> step_scores = speed_sync.search (in_data, c_speed, step, n_steps, /* seconds */ 21);
      scores.insert (scores.end(), step_scores.begin(), step_scores.end());
    }

  sort (scores.begin(), scores.end(), [] (SpeedSync::Score s_a, SpeedSync::Score s_b) { return s_a.quality > s_b.quality; });

  // we could search the N best matches, but using the best result works well in practice
  SpeedSync::Score best_s = scores[0];

  printf ("## %f %f\n", best_s.speed, best_s.quality);
  return best_s.speed;
}

static int
decode_and_report (const WavData& in_data, const string& orig_pattern)
{

  WavData wav_data;
  if (Params::detect_speed || Params::detect_speed_slow)
    {
      double speed;
      if (Params::detect_speed_slow) /* SLOW */
        {
          /* first pass:  find approximation for speed */
          speed = detect_speed (in_data, 1.0, 1.001, /* steps */ 200, /* seconds */ 15, nullptr);

          /* second pass: refine speed */
          speed = detect_speed (in_data, speed, 1.00005, /* steps */ 20,  /* seconds */ 50, nullptr);
        }
      else /* better performance, less accurate */
        {
          /* first pass:  find approximation for speed */
          speed = speed_scan (in_data);

          /* second pass: fast refine (not always perfect) */
          SpeedSync speed_sync;
          auto scores = speed_sync.search (in_data, speed, 1.00005, 20, /* seconds */ 50);
          sort (scores.begin(), scores.end(), [] (SpeedSync::Score s_a, SpeedSync::Score s_b) { return s_a.quality > s_b.quality; });
          if (!scores.empty())
            speed = scores[0].speed;
        }
      printf ("## delta %.5f %%\n", 100 * fabs (speed - Params::detect_speed_hint) / Params::detect_speed_hint);
      printf ("## speed refined %f\n", speed);

      int r = Params::mark_sample_rate * speed;
      if (r != Params::mark_sample_rate)
        wav_data = resample (in_data, Params::mark_sample_rate * speed);
      else
        wav_data = in_data;
    }
  else
    {
      wav_data = in_data;
    }
  ResultSet result_set;

  BlockDecoder block_decoder;
  block_decoder.run (wav_data, result_set);

  ClipDecoder clip_decoder;
  clip_decoder.run (wav_data, result_set);
  result_set.print();

  if (!orig_pattern.empty())
    {
      int match_count = result_set.print_match_count (orig_pattern);

      block_decoder.print_debug_sync();

      if (!match_count)
        return 1;
    }
  return 0;
}

int
get_watermark (const string& infile, const string& orig_pattern)
{
  WavData wav_data;
  Error err = wav_data.load (infile);
  if (err)
    {
      error ("audiowmark: error loading %s: %s\n", infile.c_str(), err.message());
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

/* FIXME: is this the best choice */
inline double
window_hamming (double x) /* sharp (rectangle) cutoffs at boundaries */
{
  if (fabs (x) > 1)
    return 0;

  return 0.54 + 0.46 * cos (M_PI * x);
}

void
SpeedSync::prepare_mags (const WavData& in_data, double center, double seconds)
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

  SyncFinder sync_finder;
  sync_bits = sync_finder.get_sync_bits (in_data, SyncFinder::Mode::BLOCK);

  float *in = new_array_float (sub_frame_size);
  float *out = new_array_float (sub_frame_size);

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
          fftar_float (sub_frame_size, in, out);

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

  free_array_float (in);
  free_array_float (out);
}

SpeedSync::Score
SpeedSync::compare (double relative_speed, double center)
{
  const int frames_per_block = mark_sync_frame_count() + mark_data_frame_count();
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
  return best_score;
}
