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

#include "resample.hh"
#include "wmcommon.hh"

#include <assert.h>
#include <math.h>

#include <zita-resampler/resampler.h>
#include <zita-resampler/vresampler.h>

using std::vector;

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

  WavData wav_data_out ({}, wav_data.n_channels(), rate, wav_data.bit_depth());
  vector<float>& out_ref = wav_data_out.mutable_samples();
  out_ref.resize (lrint (in.size() / wav_data.n_channels() * ratio) * wav_data.n_channels());

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
      process_resampler (resampler, in, out_ref);
      return wav_data_out;
    }

  VResampler vresampler;
  if (vresampler.setup (ratio, wav_data.n_channels(), hlen) == 0)
    {
      process_resampler (vresampler, in, out_ref);
      return wav_data_out;
    }
  error ("audiowmark: resampling from rate %d to rate %d not supported.\n", wav_data.sample_rate(), rate);
  exit (1);
}

WavData
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

template<class Resampler>
class BufferedResamplerImpl : public ResamplerImpl
{
  const int     n_channels = 0;
  const int     old_rate = 0;
  const int     new_rate = 0;
  bool          first_write = true;
  Resampler     m_resampler;

  vector<float> buffer;
public:
  BufferedResamplerImpl (int n_channels, int old_rate, int new_rate) :
    n_channels (n_channels),
    old_rate (old_rate),
    new_rate (new_rate)
  {
  }
  Resampler&
  resampler()
  {
    return m_resampler;
  }
  size_t
  skip (size_t zeros)
  {
    /* skipping a whole 1 second block should end in the same resampler state we had at the beginning */
    size_t seconds = 0;
    if (zeros >= Params::frame_size)
      seconds = (zeros - Params::frame_size) / old_rate;

    const size_t extra = new_rate * seconds;
    zeros -= old_rate * seconds;

    write_frames (vector<float> (zeros * n_channels));

    size_t out = can_read_frames() + extra;
    out -= out % Params::frame_size; /* always skip whole frames */
    read_frames (out - extra);
    return out;
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
    while (start != frames.size() / n_channels)
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

        start = frames.size() / n_channels - m_resampler.inp_count;
      }
  }
  void
  write_trailing_frames()
  {
    /* zita resampler needs k/2 samples after actual input */
    std::vector<float> samples ((m_resampler.inpsize() / 2) * n_channels);
    write_frames (samples);
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
ResamplerImpl::create (int n_channels, int old_rate, int new_rate)
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

      auto resampler = new BufferedResamplerImpl<Resampler> (n_channels, old_rate, new_rate);
      if (resampler->resampler().setup (old_rate, new_rate, n_channels, hlen) == 0)
        {
          return resampler;
        }
      else
        delete resampler;

      auto vresampler = new BufferedResamplerImpl<VResampler> (n_channels, old_rate, new_rate);
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
