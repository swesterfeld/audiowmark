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

static struct FFTPlanMap
{
  std::map<size_t, fftwf_plan> fft_plan;
  std::map<size_t, fftwf_plan> ifft_plan;

  ~FFTPlanMap()
  {
    auto free_plans = [] (auto& pmap)
      {
        for (auto it = pmap.begin(); it != pmap.end(); it++)
          {
            fftwf_destroy_plan (it->second);
          }
        pmap.clear();
      };
    free_plans (fft_plan);
    free_plans (ifft_plan);
  }
} fft_plan_map;

static std::mutex fft_planner_mutex;

FFTProcessor::FFTProcessor (size_t N)
{
  std::lock_guard<std::mutex> lg (fft_planner_mutex);

  const size_t N_2 = N + 2; /* extra space for r2c extra complex output */

  m_in  = static_cast<float *> (fftwf_malloc (sizeof (float) * N_2));
  m_out = static_cast<float *> (fftwf_malloc (sizeof (float) * N_2));

  /* plan if not done already */
  fftwf_plan& pfft = fft_plan_map.fft_plan[N];
  if (!pfft)
    pfft = fftwf_plan_dft_r2c_1d (N, m_in, (fftwf_complex *) m_out, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

  fftwf_plan& pifft = fft_plan_map.ifft_plan[N];
  if (!pifft)
    pifft = fftwf_plan_dft_c2r_1d (N, (fftwf_complex *) m_in, m_out, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

  /* store plan for size N as member variables */
  plan_fft = pfft;
  plan_ifft = pifft;

  // we could add code for saving plans here, and use patient planning
}

FFTProcessor::~FFTProcessor()
{
  fftwf_free (m_in);
  fftwf_free (m_out);
}

void
FFTProcessor::fft()
{
  fftwf_execute_dft_r2c (plan_fft, m_in, (fftwf_complex *) m_out);
}

void
FFTProcessor::ifft()
{
  fftwf_execute_dft_c2r (plan_ifft, (fftwf_complex *) m_in, m_out);
}

vector<float>
FFTProcessor::ifft (const vector<complex<float>>& in)
{
  vector<float> out ((in.size() - 1) * 2);

  /* complex<float> vector and m_out have the same layout in memory */
  std::copy (in.begin(), in.end(), reinterpret_cast<complex<float> *> (m_in));
  ifft();
  std::copy (m_out, m_out + out.size(), &out[0]);

  return out;
}

vector<complex<float>>
FFTProcessor::fft (const vector<float>& in)
{
  vector<complex<float>> out (in.size() / 2 + 1);

  /* complex<float> vector and m_out have the same layout in memory */
  std::copy (in.begin(), in.end(), m_in);
  fft();
  std::copy (m_out, m_out + out.size() * 2, reinterpret_cast<float *> (&out[0]));

  return out;
}
