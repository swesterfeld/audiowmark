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

              umag += db_from_complex (fft_out[index][u], min_db);
              dmag += db_from_complex (fft_out[index][d], min_db);
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
            umag += db_from_complex (fft_out[index][u], min_db);

          for (auto d : down)
            dmag += db_from_complex (fft_out[index][d], min_db);

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
    double            time = 0;
    vector<int>       bit_vec;
    float             decode_error = 0;
    SyncFinder::Score sync_score;
    Type              type;
    bool              speed_pattern;
  };
private:
  vector<Pattern> patterns;
  bool            speed_pattern = false;

public:
  void
  set_speed_pattern (bool sp)
  {
    speed_pattern = sp;
  }
  void
  add_pattern (double time, SyncFinder::Score sync_score, const vector<int>& bit_vec, float decode_error, Type pattern_type)
  {
    Pattern p;
    p.time = time;
    p.sync_score = sync_score;
    p.bit_vec = bit_vec;
    p.decode_error = decode_error;
    p.type = pattern_type;
    p.speed_pattern = speed_pattern;

    patterns.push_back (p);
  }
  void
  sort_by_time()
  {
    std::stable_sort (patterns.begin(), patterns.end(), [](const Pattern& p1, const Pattern& p2) {
      const int all1 = p1.type == Type::ALL;
      const int all2 = p2.type == Type::ALL;
      if (all1 != all2)
        return all1 < all2;
      else
        return p1.time < p2.time;
    });
  }
  void
  print_json (const WavData& wav_data, const std::string &json_file, const double speed)
  {
    FILE *outfile = fopen (json_file == "-" ? "/dev/stdout" : json_file.c_str(), "w");
    if (!outfile)
      {
        perror (("audiowmark: failed to open \"" + json_file + "\":").c_str());
        exit (127);
      }
    const size_t time_length = (wav_data.samples().size() / wav_data.n_channels() + wav_data.sample_rate()/2) / wav_data.sample_rate();
    fprintf (outfile, "{ \"length\": \"%ld:%02ld\",\n", time_length / 60, time_length % 60);
    fprintf (outfile, "  \"speed\": %.6f,\n", speed);
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
        if (pattern.speed_pattern)
          btype += "-SPEED";

        const int seconds = pattern.time;

        fprintf (outfile, "    { \"pos\": \"%d:%02d\", \"bits\": \"%s\", \"quality\": %.5f, \"error\": %.6f, \"type\": \"%s\" }",
                 seconds / 60, seconds % 60,
                 bit_vec_to_str (pattern.bit_vec).c_str(),
                 pattern.sync_score.quality, pattern.decode_error,
                 btype.c_str());
      }
    fprintf (outfile, " ]\n}\n");
    fclose (outfile);
  }
  void
  print()
  {
    for (const auto& pattern : patterns)
      {
        if (pattern.type == Type::ALL) /* this is the combined pattern "all" */
          {
            const char *extra = "";
            if (pattern.speed_pattern)
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
            if (pattern.speed_pattern)
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

            const double time = double (sync_score.index) / wav_data.sample_rate();
            if (!bit_vec.empty())
              result_set.add_pattern (time, sync_score, bit_vec, decode_error, ResultSet::Type::BLOCK);
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
                    result_set.add_pattern (time, score_ab, bit_vec, decode_error, ResultSet::Type::BLOCK);
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
          result_set.add_pattern (/* time */ 0.0, score_all, bit_vec, decode_error, ResultSet::Type::ALL);
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
                result_set.add_pattern (time_offset_sec, sync_score_nopad, bit_vec, decode_error, ResultSet::Type::CLIP);
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

static int
decode_and_report (const WavData& wav_data, const vector<int>& orig_bits)
{
  ResultSet result_set;
  double speed = 1.0;

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
      if (Params::detect_speed || Params::detect_speed_patient)
        speed = detect_speed (wav_data, !orig_bits.empty());
      else
        speed = Params::try_speed;

      // speeds closer to 1.0 than this usually work without stretching before decode
      if (speed < 0.9999 || speed > 1.0001)
        {
          if (Params::json_output != "-")
            printf ("speed %.6f\n", speed);
          WavData wav_data_speed = resample (wav_data, Params::mark_sample_rate * speed);

          result_set.set_speed_pattern (true);
          BlockDecoder block_decoder;
          block_decoder.run (wav_data_speed, result_set);

          ClipDecoder clip_decoder;
          clip_decoder.run (wav_data_speed, result_set);
          result_set.set_speed_pattern (false);
        }
    }

  BlockDecoder block_decoder;
  block_decoder.run (wav_data, result_set);

  ClipDecoder clip_decoder;
  clip_decoder.run (wav_data, result_set);

  result_set.sort_by_time();

  if (!Params::json_output.empty())
    result_set.print_json (wav_data, Params::json_output, speed);

  if (Params::json_output != "-")
    result_set.print();

  if (!orig_bits.empty())
    {
      int match_count = result_set.print_match_count (orig_bits);

      block_decoder.print_debug_sync();

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
get_watermark (const string& infile, const string& orig_pattern)
{
  vector<int> orig_bitvec;
  if (!orig_pattern.empty())
    {
      orig_bitvec = parse_payload (orig_pattern);
      if (orig_bitvec.empty())
        return 1;
    }

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
      return decode_and_report (wav_data, orig_bitvec);
    }
  else
    {
      return decode_and_report (resample (wav_data, Params::mark_sample_rate), orig_bitvec);
    }
}
