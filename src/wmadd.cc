#include <stdint.h>

#include <zita-resampler/resampler.h>
#include <zita-resampler/vresampler.h>

#include "wmcommon.hh"
#include "fft.hh"
#include "convcode.hh"
#include "limiter.hh"
#include "sfinputstream.hh"
#include "sfoutputstream.hh"
#include "mp3inputstream.hh"
#include "rawinputstream.hh"
#include "rawoutputstream.hh"
#include "stdoutwavoutputstream.hh"

using std::string;
using std::vector;
using std::complex;

enum class FrameMod : uint8_t {
  KEEP = 0,
  UP,
  DOWN
};

static void
prepare_frame_mod (UpDownGen& up_down_gen, int f, vector<FrameMod>& frame_mod, int data_bit)
{
  UpDownArray up, down;
  up_down_gen.get (f, up, down);
  for (auto u : up)
    frame_mod[u] = data_bit ? FrameMod::UP : FrameMod::DOWN;

  for (auto d : down)
    frame_mod[d] = data_bit ? FrameMod::DOWN : FrameMod::UP;
}

static void
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

static void
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

      for (int f = 0; f < frame_count; f++)
        {
          size_t index = data_frame_pos (f);

          prepare_frame_mod (up_down_gen, f, frame_mod[index], bitvec[f / Params::frames_per_bit]);
        }
    }
}

static void
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

static void
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

static void
init_frame_mod_vec (vector<vector<FrameMod>>& frame_mod_vec, int ab, const vector<int>& bitvec)
{
  frame_mod_vec.resize (mark_sync_frame_count() + mark_data_frame_count());

  for (auto& frame_mod : frame_mod_vec)
    frame_mod.resize (Params::max_band + 1);

  /* forward error correction */
  ConvBlockType block_type  = ab ? ConvBlockType::b : ConvBlockType::a;
  vector<int>   bitvec_fec  = randomize_bit_order (conv_encode (block_type, bitvec), /* encode */ true);

  mark_sync (frame_mod_vec, ab);
  mark_data (frame_mod_vec, bitvec_fec);
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
  int                       m_data_blocks = 0;

  FFTAnalyzer               fft_analyzer;
  WatermarkSynth            wm_synth;

  vector<int>               bitvec;
  vector<vector<FrameMod>>  pad_mod_vec;
  vector<vector<FrameMod>>  frame_mod_vec_a;
  vector<vector<FrameMod>>  frame_mod_vec_b;
public:
  WatermarkGen (int n_channels, const vector<int>& bitvec) :
    n_channels (n_channels),
    fft_analyzer (n_channels),
    wm_synth (n_channels),
    bitvec (bitvec)
  {
    init_pad_mod_vec (pad_mod_vec);
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

            // initialize a-block frame mod vector here to minimize startup latency
            assert (frame_mod_vec_a.empty());
            init_frame_mod_vec (frame_mod_vec_a, 0, bitvec);
          }
        else if (state == WATERMARK)
          {
            ab = (ab + 1) & 1; // write A|B|A|B|...
            frame_bound = mark_sync_frame_count() + mark_data_frame_count();
            m_data_blocks++;

            if (frame_mod_vec_b.empty())
              {
                // initialize b-block frame mod vector here to minimize startup latency
                init_frame_mod_vec (frame_mod_vec_b, 1, bitvec);
              }
          }
      }
    return wm_synth.run (fft_delta_spect);
  }
  int
  data_blocks() const
  {
    return m_data_blocks;
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

static ResamplerImpl *
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
  WatermarkResampler (int n_channels, int input_rate, const vector<int>& bitvec) :
    wm_gen (n_channels, bitvec)
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
  int
  data_blocks() const
  {
    return wm_gen.data_blocks();
  }
};

void
info_format (const string& label, const RawFormat& format)
{
  info ("%-13s %d Hz, %d Channels, %d Bit (%s %s-endian)\n", (label + ":").c_str(),
      format.sample_rate(), format.n_channels(), format.bit_depth(),
      format.encoding() == RawFormat::Encoding::SIGNED ? "signed" : "unsigned",
      format.endian() == RawFormat::Endian::LITTLE ? "little" : "big");
}

int
add_watermark (const string& infile, const string& outfile, const string& bits)
{
  auto bitvec = bit_str_to_vec (bits);
  if (bitvec.empty())
    {
      error ("audiowmark: cannot parse bits %s\n", bits.c_str());
      return 1;
    }
  if (bitvec.size() > Params::payload_size)
    {
      error ("audiowmark: number of bits in message '%s' larger than payload size\n", bits.c_str());
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

  /* open input stream */
  Error err;
  std::unique_ptr<AudioInputStream> in_stream = AudioInputStream::create (infile, err);
  if (err)
    {
      error ("audiowmark: error opening %s: %s\n", infile.c_str(), err.message());
      return 1;
    }

  /* open output stream */
  const int out_bit_depth = in_stream->bit_depth() > 16 ? 24 : 16;
  std::unique_ptr<AudioOutputStream> out_stream;
  out_stream = AudioOutputStream::create (outfile, in_stream->n_channels(), in_stream->sample_rate(), out_bit_depth, in_stream->n_frames(), err);
  if (err)
    {
      error ("audiowmark: error writing to %s: %s\n", outfile.c_str(), err.message());
      return 1;
    }

  /* sanity checks */
  if (in_stream->sample_rate() != out_stream->sample_rate())
    {
      error ("audiowmark: input sample rate (%d) and output sample rate (%d) don't match\n", in_stream->sample_rate(), out_stream->sample_rate());
      return 1;
    }
  if (in_stream->n_channels() != out_stream->n_channels())
    {
      error ("audiowmark: input channels (%d) and output channels (%d) don't match\n", in_stream->n_channels(), out_stream->n_channels());
      return 1;
    }

  /* write some informational messages */
  info ("Input:        %s\n", infile.c_str());
  if (Params::input_format == Format::RAW)
    info_format ("Raw Input", Params::raw_input_format);
  info ("Output:       %s\n", outfile.c_str());
  if (Params::output_format == Format::RAW)
    info_format ("Raw Output", Params::raw_output_format);
  info ("Message:      %s\n", bit_vec_to_str (bitvec).c_str());
  info ("Strength:     %.6g\n\n", Params::water_delta * 1000);

  if (in_stream->n_frames() == AudioInputStream::N_FRAMES_UNKNOWN)
    {
      info ("Time:         unknown\n");
    }
  else
    {
      int orig_seconds = in_stream->n_frames() / in_stream->sample_rate();
      info ("Time:         %d:%02d\n", orig_seconds / 60, orig_seconds % 60);
    }
  info ("Sample Rate:  %d\n", in_stream->sample_rate());
  info ("Channels:     %d\n", in_stream->n_channels());

  vector<float> samples;

  const int n_channels = in_stream->n_channels();
  AudioBuffer audio_buffer (n_channels);
  WatermarkResampler wm_resampler (n_channels, in_stream->sample_rate(), bitvec);
  if (!wm_resampler.init_ok())
    return 1;

  Limiter limiter (in_stream->sample_rate());
  limiter.set_attack (Params::limiter_attack_ms);
  limiter.set_release (Params::limiter_release_ms);
  limiter.set_ceiling (Params::limiter_ceiling);

  /* for signal to noise ratio */
  double snr_delta_power = 0;
  double snr_signal_power = 0;

  size_t total_input_frames = 0;
  size_t total_output_frames = 0;
  while (true)
    {
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

      if (Params::snr)
        {
          for (size_t i = 0; i < samples.size(); i++)
            {
              const double orig  = orig_samples[i]; // original sample
              const double delta = samples[i];      // watermark

              snr_delta_power += delta * delta;
              snr_signal_power += orig * orig;
            }
        }
      for (size_t i = 0; i < samples.size(); i++)
        samples[i] += orig_samples[i];

      samples = limiter.process (samples);

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
  if (in_stream->n_frames() != AudioInputStream::N_FRAMES_UNKNOWN)
    {
      if (total_output_frames != in_stream->n_frames())
        {
          error ("audiowmark: error: input frames (%zd) != output frames (%zd)\n", in_stream->n_frames(), total_output_frames);
          return 1;
        }
    }

  err = out_stream->close();
  if (err)
    error ("audiowmark: closing output stream failed: %s\n", err.message());

  if (Params::snr)
    info ("SNR:          %f dB\n", 10 * log10 (snr_signal_power / snr_delta_power));

  info ("Data Blocks:  %d\n", wm_resampler.data_blocks());
  return 0;
}


