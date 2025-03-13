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

#include "wavdata.hh"
#include "wmcommon.hh"
#include "wmspeed.hh"
#include "convcode.hh"
#include "shortcode.hh"
#include "syncfinder.hh"
#include "resample.hh"
#include "fft.hh"
#include "threadpool.hh"
#include "wavchunkloader.hh"

using std::string;
using std::vector;
using std::min;
using std::max;
using std::complex;

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
mix_decode (const Key& key, vector<vector<complex<float>>>& fft_out, int n_channels)
{
  vector<float> raw_bit_vec;

  const int frame_count = mark_data_frame_count();

  vector<MixEntry> mix_entries = gen_mix_entries (key);

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
              const size_t next_index = (index + n_channels) < fft_out.size() ? index + n_channels : index - n_channels;
              const size_t prev_index = (int (index) - n_channels) >= 0 ? index - n_channels : index + n_channels;

              const int u = mix_entries[b].up;
              const int d = mix_entries[b].down;

              umag += db_from_complex (fft_out[index][u], min_db);
              umag -= (db_from_complex (fft_out[prev_index][u], min_db) + db_from_complex (fft_out[next_index][u], min_db)) * 0.5;

              dmag += db_from_complex (fft_out[index][d], min_db);
              dmag -= (db_from_complex (fft_out[prev_index][d], min_db) + db_from_complex (fft_out[next_index][d], min_db)) * 0.5;
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
linear_decode (const Key& key, vector<vector<complex<float>>>& fft_out, int n_channels)
{
  UpDownGen     up_down_gen (key, Random::Stream::data_up_down);
  BitPosGen     bit_pos_gen (key);
  vector<float> raw_bit_vec;

  const int frame_count = mark_data_frame_count();

  double umag = 0, dmag = 0;
  for (int f = 0; f < frame_count; f++)
    {
      for (int ch = 0; ch < n_channels; ch++)
        {
          const size_t index = bit_pos_gen.data_frame (f) * n_channels + ch;
          const size_t next_index = (index + n_channels) < fft_out.size() ? index + n_channels : index - n_channels;
          const size_t prev_index = (int (index) - n_channels) >= 0 ? index - n_channels : index + n_channels;

          UpDownArray up, down;
          up_down_gen.get (f, up, down);

          const double min_db = -96;
          for (auto u : up)
            {
              umag += db_from_complex (fft_out[index][u], min_db);
              umag -= 0.5 * (db_from_complex (fft_out[prev_index][u], min_db) + db_from_complex (fft_out[next_index][u], min_db));
            }

          for (auto d : down)
            {
              dmag += db_from_complex (fft_out[index][d], min_db);
              dmag -= 0.5 * (db_from_complex (fft_out[prev_index][d], min_db) + db_from_complex (fft_out[next_index][d], min_db));
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
mix_or_linear_decode (const Key& key, vector<vector<complex<float>>>& fft_out, int n_channels)
{
  if (Params::mix)
    return mix_decode (key, fft_out, n_channels);
  else
    return linear_decode (key, fft_out, n_channels);
}

class ResultSet
{
public:
  enum class Type { BLOCK, CLIP, ALL };
  struct Pattern
  {
    Key               key;
    double            time = 0;
    vector<int>       bit_vec;
    float             decode_error = 0;
    SyncFinder::Score sync_score;
    Type              type;
    double            speed = 0;

    bool
    approx_match (const Pattern& p) const
    {
      const double time_delta = Params::frame_size / double (Params::mark_sample_rate);
      const double speed_delta = 0.01;

      return key == p.key &&
            (fabs (time - p.time) < time_delta || (type == Type::ALL)) &&
             bit_vec == p.bit_vec &&
             sync_score.block_type == p.sync_score.block_type &&
             type == p.type &&
             fabs (speed - p.speed) < speed_delta;
    }
  };
private:
  std::mutex      pattern_mutex;
  vector<Pattern> patterns;
  std::string     debug_sync;

public:
  void
  add_pattern (const Key& key, double time, SyncFinder::Score sync_score, const vector<int>& bit_vec, float decode_error, Type pattern_type, double speed)
  {
    /* add_pattern can be called by any thread (safe to use from ThreadPool jobs) */
    std::lock_guard<std::mutex> lg (pattern_mutex);

    Pattern p;
    p.key = key;
    p.time = time;
    p.sync_score = sync_score;
    p.bit_vec = bit_vec;
    p.decode_error = decode_error;
    p.type = pattern_type;
    p.speed = speed;

    patterns.push_back (p);
  }
  void
  apply_time_offset (double time_offset)
  {
    for (auto& p : patterns)
      p.time += time_offset;
  }
  void
  sort()
  {
    std::sort (patterns.begin(), patterns.end(), [](const Pattern& p1, const Pattern& p2) {
      const int all1 = p1.type == Type::ALL;
      const int all2 = p2.type == Type::ALL;
      const auto p1bits = bit_vec_to_str (p1.bit_vec);
      const auto p2bits = bit_vec_to_str (p2.bit_vec);

      auto ab = [] (const Pattern& pattern) {
        switch (pattern.sync_score.block_type) {
          case ConvBlockType::a:  return 0;
          case ConvBlockType::b:  return 1;
          case ConvBlockType::ab: return 2;
        };
        return 99; // should not happen
      };

      if (p1.key.name() != p2.key.name())
        return p1.key.name() < p2.key.name();
      else if (all1 != all2)
        return all1 < all2;
      else if (p1.time != p2.time)
        return p1.time < p2.time;
      else if (ab (p1) != ab (p2))
        return ab (p1) < ab (p2);
      else
        {
          return p1bits < p2bits;
        }
    });
  }
  void
  merge (ResultSet& other)
  {
    for (const auto& p : other.patterns)
      {
        bool merge = true;
        for (auto& my_p : patterns)
          {
            if (my_p.approx_match (p))
              merge = false;
          }
        if (merge)
          patterns.push_back (p);
      }

    /* only keep track of debug sync information for the first chunk */
    if (debug_sync.empty())
      debug_sync = other.debug_sync;
  }
  string
  json_escape (const string& s)
  {
    string result;
    for (unsigned ch : s)
      {
        if (ch == '"' || ch == '\\')
          {
            result += '\\';
            result += ch;
          }
        else if (ch < 32)
          {
            result += string_printf ("\\u%04x", ch);
          }
        else
          {
            result += ch;
          }
      }
    return result;
  }
  void
  print_json (size_t time_length, const std::string &json_file)
  {
    FILE *outfile = fopen (json_file == "-" ? "/dev/stdout" : json_file.c_str(), "w");
    if (!outfile)
      {
        perror (("audiowmark: failed to open \"" + json_file + "\":").c_str());
        exit (127);
      }
    fprintf (outfile, "{ \"length\": \"%ld:%02ld\",\n", time_length / 60, time_length % 60);
    fprintf (outfile, "  \"matches\": [\n");
    int nth = 0;
    for (const auto& pattern : patterns)
      {
        if (nth++ != 0)
          fprintf (outfile, ",\n");

        std::string btype;
        switch (pattern.sync_score.block_type)
          {
          case ConvBlockType::a:        btype = "A";    break;
          case ConvBlockType::b:        btype = "B";    break;
          case ConvBlockType::ab:       btype = "AB";   break;
          }
        if (pattern.type == Type::ALL)
          btype = "ALL";
        if (pattern.type == Type::CLIP)
          btype = "CLIP-" + btype;
        if (pattern.speed != 1)
          btype += "-SPEED";

        const int seconds = pattern.time;

        fprintf (outfile, "    { \"key\": \"%s\", \"pos\": \"%d:%02d\", \"bits\": \"%s\", \"quality\": %.5f, \"error\": %.6f, \"type\": \"%s\", \"speed\": %.6f }",
                 json_escape (pattern.key.name()).c_str(),
                 seconds / 60, seconds % 60,
                 bit_vec_to_str (pattern.bit_vec).c_str(),
                 pattern.sync_score.quality, pattern.decode_error,
                 btype.c_str(),
                 pattern.speed);
      }
    fprintf (outfile, " ]\n}\n");
    fclose (outfile);
  }
  void
  print()
  {
    string last_key_name;

    for (const auto& pattern : patterns)
      {
        if (pattern.key.name() != last_key_name)
          {
            printf ("key %s\n", pattern.key.name().c_str());
            last_key_name = pattern.key.name();

            // currently we assume that speed detection returns one best speed for each key
            for (auto p : patterns)
              if (p.key.name() == pattern.key.name() && p.speed != 1)
                {
                  printf ("speed %.6f\n", p.speed);
                  break;
                }
          }
        if (pattern.type == Type::ALL) /* this is the combined pattern "all" */
          {
            const char *extra = "";
            if (pattern.speed != 1)
              extra = " SPEED";

            printf ("pattern   all %s %.3f %.3f%s\n", bit_vec_to_str (pattern.bit_vec).c_str(),
                                                      pattern.sync_score.quality, pattern.decode_error,
                                                      extra);
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
            if (pattern.speed != 1)
              block_str += "-SPEED";

            const int seconds = pattern.time;
            printf ("pattern %2d:%02d %s %.3f %.3f %s\n", seconds / 60, seconds % 60, bit_vec_to_str (pattern.bit_vec).c_str(),
                      pattern.sync_score.quality, pattern.decode_error, block_str.c_str());
          }
      }
  }
  int
  print_match_count (const vector<int>& orig_bits)
  {
    int match_count = 0;

    for (auto p : patterns)
      {
        if (p.bit_vec == orig_bits)
          match_count++;
      }
    printf ("match_count %d %zd\n", match_count, patterns.size());
    return match_count;
  }
  void
  set_debug_sync (const std::string& ds)
  {
    debug_sync = ds;
  }
  void
  print_debug_sync()
  {
    printf ("%s", debug_sync.c_str());
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
  const double speed = 0;
  vector<SyncFinder::KeyResult> key_results; // stored here for sync debugging
public:
  BlockDecoder (double speed) :
    speed (speed)
  {
  }
  void
  run (const vector<Key>& key_list, const WavData& wav_data, ResultSet& result_set)
  {
    ThreadPool thread_pool;
    SyncFinder sync_finder;
    FFTAnalyzer fft_analyzer (wav_data.n_channels());
    key_results = sync_finder.search (key_list, wav_data, SyncFinder::Mode::BLOCK);

    for (const auto& key_result : key_results)
      {
        const Key&  key = key_result.key;

        struct PatternRawBits {
          size_t        index;
          double        quality;
          vector<float> raw_bit_vec;
          ConvBlockType block_type;
        };
        vector<PatternRawBits> pattern_raw_vec;
        for (auto sync_score : key_result.sync_scores)
          {
            const size_t count = mark_sync_frame_count() + mark_data_frame_count();
            const size_t index = sync_score.index;

            auto fft_range_out = fft_analyzer.fft_range (wav_data.samples(), index, count);
            if (fft_range_out.size())
              {
                /* ---- retrieve bits from watermark ---- */
                vector<float> raw_bit_vec = mix_or_linear_decode (key, fft_range_out, wav_data.n_channels());
                assert (raw_bit_vec.size() == code_size (ConvBlockType::a, Params::payload_size));

                raw_bit_vec = randomize_bit_order (key, raw_bit_vec, /* encode */ false);

                PatternRawBits raw_bits;
                raw_bits.index = index;
                raw_bits.quality = sync_score.quality;
                raw_bits.raw_bit_vec = raw_bit_vec;
                raw_bits.block_type = sync_score.block_type;
                pattern_raw_vec.push_back (raw_bits);

                /* ---- deal with this pattern ---- */
                const double time = double (sync_score.index) / wav_data.sample_rate();
                thread_pool.add_job ([this, key, sync_score, raw_bit_vec, time, &result_set]()
                  {
                    float decode_error = 0;
                    vector<int> bit_vec = code_decode_soft (sync_score.block_type, normalize_soft_bits (raw_bit_vec), &decode_error);

                    if (!bit_vec.empty())
                      result_set.add_pattern (key, time, sync_score, bit_vec, decode_error, ResultSet::Type::BLOCK, speed);
                  });
              }
          }
        /* AB pattern: try to find an A block followed by a B block with the right distance (sync + data frame count) */
        for (size_t i = 0; i < pattern_raw_vec.size(); i++)
          {
            if (pattern_raw_vec[i].block_type == ConvBlockType::b)
              {
                int best_j = -1;
                int best_abs_dist = Params::frame_size / 2;
                for (size_t j = 0; j < i; j++)
                  {
                    if (pattern_raw_vec[j].block_type == ConvBlockType::a)
                      {
                        size_t count = mark_sync_frame_count() + mark_data_frame_count();
                        int    abs_dist = std::abs (int (pattern_raw_vec[i].index - pattern_raw_vec[j].index) - int (count * Params::frame_size));

                        if (abs_dist < best_abs_dist)
                          {
                            best_j = j;
                            best_abs_dist = abs_dist;
                          }
                      }
                  }
                /* join A and B block -> AB block */
                if (best_j >= 0)
                  {
                    const auto& a_pattern = pattern_raw_vec[best_j];
                    const auto& b_pattern = pattern_raw_vec[i];

                    vector<float> ab_bits (a_pattern.raw_bit_vec.size() * 2);
                    for (size_t k = 0; k <  a_pattern.raw_bit_vec.size(); k++)
                      {
                        ab_bits[k * 2]     = a_pattern.raw_bit_vec[k];
                        ab_bits[k * 2 + 1] = b_pattern.raw_bit_vec[k];
                      }

                    const double time = double (b_pattern.index) / wav_data.sample_rate();
                    thread_pool.add_job ([this, key, a_pattern, b_pattern, ab_bits, time, &result_set]()
                      {
                        float decode_error = 0;
                        vector<int> bit_vec = code_decode_soft (ConvBlockType::ab, normalize_soft_bits (ab_bits), &decode_error);

                        if (!bit_vec.empty())
                          {
                            SyncFinder::Score score_ab  { 0, 0, ConvBlockType::ab };
                            score_ab.index = b_pattern.index;
                            score_ab.quality = (a_pattern.quality + b_pattern.quality) / 2;
                            result_set.add_pattern (key, time, score_ab, bit_vec, decode_error, ResultSet::Type::BLOCK, speed);
                          }
                      });
                  }
              }
          }

        /* all pattern: try to detect consecutive blocks with the right distance (sync + data frame count) */
        vector<size_t> best_all_blocks;
        for (size_t i = 0; i < pattern_raw_vec.size(); i++)
          {
            vector<size_t> all_blocks;
            all_blocks.push_back (i);
            for (size_t block_idx = 1; block_idx < pattern_raw_vec.size(); block_idx++)
              {
                size_t count = mark_sync_frame_count() + mark_data_frame_count();
                size_t expect_start = pattern_raw_vec[i].index + block_idx * (count * Params::frame_size);
                int    best_j = -1;
                int    best_abs_dist = Params::frame_size / 2;
                for (size_t j = i; j < pattern_raw_vec.size(); j++)
                  {
                    int abs_dist = std::abs (int (expect_start) - int (pattern_raw_vec[j].index));
                    if (abs_dist < best_abs_dist)
                      {
                        best_j = j;
                        best_abs_dist = abs_dist;
                      }
                  }

                if (best_j >= 0)
                  all_blocks.push_back (best_j);
              }
            if (all_blocks.size() == best_all_blocks.size())
              {
                /* for two all patterns which have the same amount of blocks:
                 * prefer pattern with higher sync score sum
                 */
                auto sync_sum = [&] (auto blocks)
                  {
                    float sum = 0;
                    for (auto block_idx : blocks)
                      sum += pattern_raw_vec[block_idx].quality;
                    return sum;
                  };
                if (sync_sum (all_blocks) > sync_sum (best_all_blocks))
                    best_all_blocks = all_blocks;
              }
            if (all_blocks.size() > best_all_blocks.size())
              {
                /* prefer all patterns with higher number of blocks */
                best_all_blocks = all_blocks;
              }
          }
        /* all pattern: average the A / B bits of the consecutive blocks for an "all" pattern */
        if (best_all_blocks.size() > 1)
          {
            vector<float> raw_bit_vec_all (code_size (ConvBlockType::ab, Params::payload_size));
            vector<int>   raw_bit_vec_norm (2);

            SyncFinder::Score score_all { 0, 0 };

            for (auto block_index : best_all_blocks)
              {
                const auto& pattern = pattern_raw_vec[block_index];
                /* ---- update "all" pattern ---- */
                score_all.quality += pattern.quality;

                int ab = pattern.block_type == ConvBlockType::b ? 1 : 0;
                for (size_t i = 0; i < pattern.raw_bit_vec.size(); i++)
                  {
                    raw_bit_vec_all[i * 2 + ab] += pattern.raw_bit_vec[i];
                  }
                raw_bit_vec_norm[ab]++;
              }

            for (size_t i = 0; i < raw_bit_vec_all.size(); i += 2)
              {
                raw_bit_vec_all[i]     /= max (raw_bit_vec_norm[0], 1); /* normalize A soft bits with number of A blocks */
                raw_bit_vec_all[i + 1] /= max (raw_bit_vec_norm[1], 1); /* normalize B soft bits with number of B blocks */
              }
            score_all.quality /= raw_bit_vec_norm[0] + raw_bit_vec_norm[1];

            vector<float> soft_bit_vec = normalize_soft_bits (raw_bit_vec_all);

            thread_pool.add_job ([this, key, score_all, soft_bit_vec, &result_set]()
              {
                float decode_error = 0;
                vector<int> bit_vec = code_decode_soft (ConvBlockType::ab, soft_bit_vec, &decode_error);

                if (!bit_vec.empty())
                  result_set.add_pattern (key, /* time */ 0.0, score_all, bit_vec, decode_error, ResultSet::Type::ALL, speed);
              });
          }
      }
    thread_pool.wait_all();

    debug_sync_frame_count = frame_count (wav_data);
  }
  std::string
  debug_sync()
  {
    /* this is really only useful for debugging, and should be used with exactly one key */
    if (key_results.size() != 1)
      return "";

    const auto& sync_scores = key_results[0].sync_scores;

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
    return string_printf ("sync_match %d %zd\n", sync_match, sync_scores.size());
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
  const double speed = 0;

  void
  run_padded (const vector<Key>& key_list, const WavData& wav_data, ResultSet& result_set, double time_offset_sec)
  {
    SyncFinder                    sync_finder;
    vector<SyncFinder::KeyResult> key_results = sync_finder.search (key_list, wav_data, SyncFinder::Mode::CLIP);
    FFTAnalyzer                   fft_analyzer (wav_data.n_channels());
    ThreadPool                    thread_pool;

    for (const auto& key_result : key_results)
      {
        const Key& key = key_result.key;
        for (const auto& sync_score : key_result.sync_scores)
          {
            const size_t count = mark_sync_frame_count() + mark_data_frame_count();
            const size_t index = sync_score.index;
            auto fft_range_out1 = fft_analyzer.fft_range (wav_data.samples(), index, count);
            auto fft_range_out2 = fft_analyzer.fft_range (wav_data.samples(), index + count * Params::frame_size, count);
            if (fft_range_out1.size() && fft_range_out2.size())
              {
                const auto raw_bit_vec1 = randomize_bit_order (key, mix_or_linear_decode (key, fft_range_out1, wav_data.n_channels()), /* encode */ false);
                const auto raw_bit_vec2 = randomize_bit_order (key, mix_or_linear_decode (key, fft_range_out2, wav_data.n_channels()), /* encode */ false);
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

                SyncFinder::Score sync_score_nopad = sync_score;
                sync_score_nopad.index = time_offset_sec * wav_data.sample_rate();

                thread_pool.add_job ([this, key, raw_bit_vec, sync_score_nopad, time_offset_sec, &result_set]()
                  {
                    float decode_error = 0;
                    vector<int> bit_vec = code_decode_soft (ConvBlockType::ab, normalize_soft_bits (raw_bit_vec), &decode_error);

                    if (!bit_vec.empty())
                      result_set.add_pattern (key, time_offset_sec, sync_score_nopad, bit_vec, decode_error, ResultSet::Type::CLIP, speed);
                  });
              }
          }
      }
    thread_pool.wait_all();
  }
  enum class Pos { START, END };
  void
  run_block (const vector<Key>& key_list, const WavData& wav_data, ResultSet& result_set, Pos pos)
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
    run_padded (key_list, l_wav_data, result_set, time_offset);
   }
public:
  ClipDecoder(double speed) :
    frames_per_block (mark_sync_frame_count() + mark_data_frame_count()),
    speed (speed)
  {
  }
  void
  run (const vector<Key>& key_list, const WavData& wav_data, ResultSet& result_set)
  {
    const int wav_frames = wav_data.n_values() / (Params::frame_size * wav_data.n_channels());
    if (wav_frames < frames_per_block * 3.1) /* clip decoder is only used for small wavs */
      {
        run_block (key_list, wav_data, result_set, Pos::START);
        run_block (key_list, wav_data, result_set, Pos::END);
      }
  }
};

static void
decode (ResultSet& result_set, const vector<Key>& key_list, const WavData& wav_data, const vector<int>& orig_bits, bool first_chunk)
{
  /*
   * The strategy for integrating speed detection into decoding is this:
   *  - we always (unconditionally) try to decode the  watermark on the original wav data
   *  - if detected speed is somewhat different than 1.0, we also try to decode stretched data
   *  - we report all normal and speed results we get
   *
   * The reason to do it this way is that the detected speed may be wrong (on short clips)
   * and we don't want to loose a successful clip decoder match in this case.
   */
  if (Params::detect_speed || Params::detect_speed_patient || Params::try_speed > 0)
    {
      vector<DetectSpeedResult> speed_results;
      if (Params::detect_speed || Params::detect_speed_patient)
        speed_results = detect_speed (key_list, wav_data, !orig_bits.empty());
      else
        {
          for (const auto& key : key_list)
            {
              DetectSpeedResult speed_result;
              speed_result.key   = key;
              speed_result.speed = Params::try_speed;
              speed_results.push_back (speed_result);
            }
        }

      for (const auto& speed_result : speed_results)
        {
          WavData wav_data_speed = resample_ratio (wav_data, speed_result.speed, Params::mark_sample_rate * speed_result.speed);

          BlockDecoder block_decoder (speed_result.speed);
          block_decoder.run ({ speed_result.key }, wav_data_speed, result_set);

          if (first_chunk)
            {
              ClipDecoder clip_decoder (speed_result.speed);
              clip_decoder.run ({ speed_result.key }, wav_data_speed, result_set);
            }
        }
    }

  BlockDecoder block_decoder (1);
  block_decoder.run (key_list, wav_data, result_set);

  if (first_chunk)
    {
      ClipDecoder clip_decoder (1);
      clip_decoder.run (key_list, wav_data, result_set);
    }

  result_set.set_debug_sync (block_decoder.debug_sync());
}

int
report (ResultSet& result_set, size_t time_length, const vector<int>& orig_bits)
{
  if (!Params::json_output.empty())
    result_set.print_json (time_length, Params::json_output);

  if (Params::json_output != "-")
    result_set.print();

  if (!orig_bits.empty())
    {
      int match_count = result_set.print_match_count (orig_bits);

      result_set.print_debug_sync();

      if (Params::expect_matches >= 0)
        {
          printf ("expect_matches %d\n", Params::expect_matches);
          if (match_count != Params::expect_matches)
            return 1;
        }
      else
        {
          if (!match_count)
            return 1;
        }
    }
  return 0;
}

int
get_watermark (const vector<Key>& key_list, const string& infile, const string& orig_pattern)
{
  ResultSet result_set;

  vector<int> orig_bitvec;
  if (!orig_pattern.empty())
    {
      orig_bitvec = parse_payload (orig_pattern);
      if (orig_bitvec.empty())
        return 1;
    }

  bool first_chunk = true;
  WavChunkLoader wav_chunk_loader (infile);
  while (!wav_chunk_loader.done())
    {
      Error err = wav_chunk_loader.load_next_chunk();
      if (err)
        {
          error ("audiowmark: error loading %s: %s\n", infile.c_str(), err.message());
          return 1;
        }

      if (!wav_chunk_loader.done())
        {
          const WavData& wav_data = wav_chunk_loader.wav_data();
          assert (wav_data.sample_rate() == Params::mark_sample_rate);

          ResultSet chunk_result_set;

          decode (chunk_result_set, key_list, wav_data, orig_bitvec, first_chunk);
          chunk_result_set.apply_time_offset (wav_chunk_loader.time_offset());

          result_set.merge (chunk_result_set);
          first_chunk = false;
        }
    }
  result_set.sort();

  size_t time_length = lrint (wav_chunk_loader.length());
  return report (result_set, time_length, orig_bitvec);
}
