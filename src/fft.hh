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

#ifndef AUDIOWMARK_FFT_HH
#define AUDIOWMARK_FFT_HH

#include <complex>
#include <vector>
#include <fftw3.h>

class FFTProcessor
{
  fftwf_plan plan_fft;
  fftwf_plan plan_ifft;
  float *m_in = nullptr;
  float *m_out = nullptr;
public:
  FFTProcessor (size_t N);
  ~FFTProcessor();

  /* low level (fast) */
  void   fft();
  void   ifft();
  float *in()  { return m_in; }
  float *out() { return m_out; };

  /* high level (convenient) */
  std::vector<std::complex<float>> fft (const std::vector<float>& in);
  std::vector<float>               ifft (const std::vector<std::complex<float>>& in);
};

#endif /* AUDIOWMARK_FFT_HH */
