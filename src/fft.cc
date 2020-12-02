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

#include "fft.hh"

#include <fftw3.h>

#include <map>
#include <mutex>

using std::vector;
using std::complex;
using std::map;

float *
new_array_float (size_t N)
{
  const size_t N_2 = N + 2; /* extra space for r2c extra complex output */

  return (float *) fftwf_malloc (sizeof (float) * N_2);
}

void
free_array_float (float *f)
{
  fftwf_free (f);
}

static std::mutex fft_planner_mutex;

FFTProcessor::FFTProcessor (size_t N)
{
  std::lock_guard<std::mutex> lg (fft_planner_mutex);

  float *plan_in = new_array_float (N);
  float *plan_out = new_array_float (N);

  plan_fft = fftwf_plan_dft_r2c_1d (N, plan_in, (fftwf_complex *) plan_out, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  plan_ifft = fftwf_plan_dft_c2r_1d (N, (fftwf_complex *) plan_in, plan_out, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

  // we add code for saving plans here, and use patient planning

  free_array_float (plan_out);
  free_array_float (plan_in);
}

void
FFTProcessor::fft (float *in, float *out)
{
  fftwf_execute_dft_r2c (plan_fft, in, (fftwf_complex *) out);
}

void
FFTProcessor::ifft (float *in, float *out)
{
  fftwf_execute_dft_c2r (plan_ifft, (fftwf_complex *)in, out);
}

vector<float>
FFTProcessor::ifft (const vector<complex<float>>& in)
{
  vector<float> out ((in.size() - 1) * 2);

  /* ensure memory is SSE-aligned (or other vectorized stuff) */
  float *ifft_in = new_array_float (out.size());
  float *ifft_out = new_array_float (out.size());

  /* complex<float> vector and fft_out have the same layout in memory */
  std::copy (in.begin(), in.end(), reinterpret_cast<complex<float> *> (ifft_in));
  ifft (ifft_in, ifft_out);

  std::copy (ifft_out, ifft_out + out.size(), &out[0]);

  free_array_float (ifft_out);
  free_array_float (ifft_in);

  return out;
}

vector<complex<float>>
FFTProcessor::fft (const vector<float>& in)
{
  vector<complex<float>> out (in.size() / 2 + 1);

  /* ensure memory is SSE-aligned (or other vectorized stuff) */
  float *fft_in = new_array_float (in.size());
  float *fft_out = new_array_float (in.size());

  std::copy (in.begin(), in.end(), fft_in);
  fft (fft_in, fft_out);

  /* complex<float> vector and fft_out have the same layout in memory */
  std::copy (fft_out, fft_out + out.size() * 2, reinterpret_cast<float *> (&out[0]));

  free_array_float (fft_out);
  free_array_float (fft_in);

  return out;
}
