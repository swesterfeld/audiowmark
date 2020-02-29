#include <string>
#include <algorithm>

#include <zita-resampler/resampler.h>
#include <zita-resampler/vresampler.h>

#include "wavdata.hh"
#include "wmcommon.hh"
#include "convcode.hh"

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

static int
frame_count (const WavData& wav_data)
{
  return wav_data.n_values() / wav_data.n_channels() / Params::frame_size;
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

static double
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
public:
  enum class BlockLength { NORMAL, LONG };

private:
  struct FrameBit
  {
    int frame;
    vector<int> up;
    vector<int> down;
  };
  vector<vector<FrameBit>> sync_bits;

  void
  init_up_down (const WavData& wav_data, BlockLength block_length)
  {
    sync_bits.clear();

    // "long" blocks consist of two "normal" blocks, which means
    //   the sync bits pattern is repeated after the end of the first block
    const int first_block_end = mark_sync_frame_count() + mark_data_frame_count();
    const int block_count = block_length == BlockLength::LONG ? 2 : 1;
    size_t n_bands = Params::max_band - Params::min_band + 1;
    for (int block = 0; block < block_count; block++)
      {
        UpDownGen up_down_gen (Random::Stream::sync_up_down);
        for (int bit = 0; bit < Params::sync_bits; bit++)
          {
            vector<FrameBit> frame_bits;
            for (int f = 0; f < Params::sync_frames_per_bit; f++)
              {
                UpDownArray frame_up, frame_down;
                up_down_gen.get (f + bit * Params::sync_frames_per_bit, frame_up, frame_down);

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
            std::sort (frame_bits.begin(), frame_bits.end(), [] (FrameBit& f1, FrameBit& f2) { return f1.frame < f2.frame; });
            sync_bits.push_back (frame_bits);
          }
      }
  }
  double
  sync_decode (const WavData& wav_data, const size_t start_frame, const vector<float>& fft_out_db, ConvBlockType *block_type, int *len)
  {
    double sync_quality = 0;

    size_t n_bands = Params::max_band - Params::min_band + 1;
    int lbits = 0;
    for (size_t bit = 0; bit < sync_bits.size(); bit++)
      {
        const vector<FrameBit>& frame_bits = sync_bits[bit];
        float umag = 0, dmag = 0;

        for (const auto& frame_bit : frame_bits)
          {
            int index = ((start_frame + frame_bit.frame) * wav_data.n_channels()) * n_bands;
            if (fft_out_db[index] > -90)
              {
                for (size_t i = 0; i < frame_bit.up.size(); i++)
                  {
                    umag += fft_out_db[index + frame_bit.up[i]];
                    dmag += fft_out_db[index + frame_bit.down[i]];
                  }
                lbits++;
              }
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

        const int expect_data_bit = bit & 1; /* expect 010101 */
        const double q = expect_data_bit ? raw_bit : -raw_bit;
        sync_quality += q;
      }
    sync_quality /= sync_bits.size();
    sync_quality = normalize_sync_quality (sync_quality);

    *len = lbits;
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
    int           len = 0;
  };
  vector<Score>
  search (const WavData& wav_data, BlockLength block_length)
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
    init_up_down (wav_data, block_length);

    vector<float> fft_db;

    // compute multiple time-shifted fft vectors
    size_t n_bands = Params::max_band - Params::min_band + 1;
    int total_frame_count = mark_sync_frame_count() + mark_data_frame_count();
    const int first_block_end = total_frame_count;
    if (block_length == BlockLength::LONG)
      total_frame_count *= 2;
    for (size_t sync_shift = 0; sync_shift < Params::frame_size; sync_shift += Params::sync_search_step)
      {
        sync_fft (wav_data, sync_shift, frame_count (wav_data) - 1, fft_db, /* want all frames */ {});
        for (int start_frame = 0; start_frame < frame_count (wav_data); start_frame++)
          {
            const size_t sync_index = start_frame * Params::frame_size + sync_shift;
            if ((start_frame + total_frame_count) * wav_data.n_channels() * n_bands < fft_db.size())
              {
                ConvBlockType block_type;
                int len;
                double quality = sync_decode (wav_data, start_frame, fft_db, &block_type, &len);
                // printf ("%zd %f\n", sync_index, quality);
                sync_scores.emplace_back (Score { sync_index, quality, block_type, len });
              }
          }
      }
    sort (sync_scores.begin(), sync_scores.end(), [] (const Score& a, const Score &b) { return a.index < b.index; });

    vector<int> want_frames (total_frame_count);
    for (size_t f = 0; f < mark_sync_frame_count(); f++)
      {
        want_frames[sync_frame_pos (f)] = 1;
        if (block_length == BlockLength::LONG)
          want_frames[first_block_end + sync_frame_pos (f)] = 1;
      }

    /* for strength 8 and above:
     *   -> more false positive candidates are rejected, so we can use a lower threshold
     *
     * for strength 7 and below:
     *   -> we need a higher threshold, because otherwise watermark detection takes too long
     */
    const double strength = Params::water_delta * 1000;
    const double sync_threshold1 = strength > 7.5 ? 0.4 : 0.5;

    /* n-best-search */
    struct NBest {
      size_t i;
      int len;
    };
    vector<NBest> n_best;
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

            if (sync_scores[i].quality > q_last && sync_scores[i].quality > q_next)
              {
                n_best.push_back ({ i, sync_scores[i].len });
              }
          }
      }
#if 0
    std::sort (n_best.begin(), n_best.end(), [](NBest& nb1, NBest& nb2) { return nb1.len > nb2.len; });
    if (n_best.size() > 5)
      n_best.resize (5);
#endif
    for (auto n : n_best)
      {
        size_t i = n.i;
        //printf ("%d %zd %f %d\n", int (block_length), sync_scores[i].index, sync_scores[i].quality, sync_scores[i].len);
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
                    sync_fft (wav_data, fine_index, total_frame_count, fft_db, want_frames);
                    if (fft_db.size())
                      {
                        ConvBlockType block_type;
                        int           len;
                        double        q = sync_decode (wav_data, 0, fft_db, &block_type, &len);

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
            bool zframe = true;
            for (int i = 0; i < Params::frame_size * wav_data.n_channels(); i++)
              {
                if (samples[i + (index + f * Params::frame_size) * wav_data.n_channels()] != 0)
                  {
                    zframe = false;
                    break;
                  }
              }
            if (zframe)
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

static int
decode_and_report (const WavData& wav_data, const string& orig_pattern)
{
  int match_count = 0, total_count = 0, sync_match = 0;

  SyncFinder                sync_finder;
  vector<SyncFinder::Score> sync_scores; // = sync_finder.search (wav_data, SyncFinder::BlockLength::NORMAL);

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
  // L-blocks
  SyncFinder                sync_finder_l;
  vector<SyncFinder::Score> sync_scores_l = sync_finder.search (wav_data, SyncFinder::BlockLength::LONG);

  for (auto sync_score : sync_scores_l)
    {
      string block_str;
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
      const size_t count = mark_sync_frame_count() + mark_data_frame_count();
      const size_t index = sync_score.index;
      auto fft_range_out1 = fft_analyzer.fft_range (wav_data.samples(), index, count);
      auto fft_range_out2 = fft_analyzer.fft_range (wav_data.samples(), index + count * Params::frame_size, count);
      if (fft_range_out1.size() && fft_range_out2.size())
        {
          // FIXME: doesn't do linear decode
          const auto raw_bit_vec1 = randomize_bit_order (mix_decode (fft_range_out1, wav_data.n_channels()), /* encode */ false);
          const auto raw_bit_vec2 = randomize_bit_order (mix_decode (fft_range_out2, wav_data.n_channels()), /* encode */ false);
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
          vector<int> bit_vec = conv_decode_soft (ConvBlockType::ab, normalize_soft_bits (raw_bit_vec), &decode_error);
          printf ("LLL ");
          report_pattern (sync_score, bit_vec, decode_error);
        }
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
  auto ext_samples = wav_data.samples();
  ext_samples.insert (ext_samples.begin(), wav_data.sample_rate() * wav_data.n_channels() * 100, 0);
  ext_samples.insert (ext_samples.end(),   wav_data.sample_rate() * wav_data.n_channels() * 100, 0);
  wav_data.set_samples (ext_samples);
  if (wav_data.sample_rate() == Params::mark_sample_rate)
    {
      return decode_and_report (wav_data, orig_pattern);
    }
  else
    {
      return decode_and_report (resample (wav_data, Params::mark_sample_rate), orig_pattern);
    }
}
