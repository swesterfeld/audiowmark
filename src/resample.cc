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


