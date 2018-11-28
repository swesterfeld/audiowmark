#ifndef AUDIOWMARK_FFT_HH
#define AUDIOWMARK_FFT_HH

#include <complex>
#include <vector>

std::vector<std::complex<float>> fft (const std::vector<float>& in);
std::vector<float>               ifft (const std::vector<std::complex<float>>& in);

#endif /* AUDIOWMARK_FFT_HH */
