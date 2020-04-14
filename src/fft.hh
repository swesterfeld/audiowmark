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

/* high level api */
std::vector<std::complex<float>> fft (const std::vector<float>& in);
std::vector<float>               ifft (const std::vector<std::complex<float>>& in);

/* more efficient: low level api */
void   fftar_float (size_t N, float *in, float *out);
float *new_array_float (size_t N);
void   free_array_float (float *f);


#endif /* AUDIOWMARK_FFT_HH */
