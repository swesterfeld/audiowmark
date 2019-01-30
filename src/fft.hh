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
