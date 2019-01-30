#include "fft.hh"

#include <fftw3.h>

#include <map>

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

void
fftar_float (size_t N, float *in, float *out)
{
  static map<int, fftwf_plan> plan_for_size;

  fftwf_plan& plan = plan_for_size[N];
  if (!plan)
    {
      float *plan_in = new_array_float (N);
      float *plan_out = new_array_float (N);
      plan = fftwf_plan_dft_r2c_1d (N, plan_in, (fftwf_complex *) plan_out, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

      // we add code for saving plans here, and use patient planning
    }
  fftwf_execute_dft_r2c (plan, in, (fftwf_complex *) out);
}

void
fftsr_float (size_t N, float *in, float *out)
{
  static map<int, fftwf_plan> plan_for_size;

  fftwf_plan& plan = plan_for_size[N];
  if (!plan)
    {
      float *plan_in = new_array_float (N);
      float *plan_out = new_array_float (N);
      plan = fftwf_plan_dft_c2r_1d (N, (fftwf_complex *) plan_in, plan_out, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

      // we add code for saving plans here, and use patient planning
    }
  fftwf_execute_dft_c2r (plan, (fftwf_complex *)in, out);
}

vector<complex<float>>
fft (const vector<float>& in)
{
  vector<complex<float>> out (in.size() / 2 + 1);

  /* ensure memory is SSE-aligned (or other vectorized stuff) */
  float *fft_in = new_array_float (in.size());
  float *fft_out = new_array_float (in.size());

  std::copy (in.begin(), in.end(), fft_in);
  fftar_float (in.size(), fft_in, fft_out);

  /* complex<float> vector and fft_out have the same layout in memory */
  std::copy (fft_out, fft_out + out.size() * 2, reinterpret_cast<float *> (&out[0]));

  free_array_float (fft_out);
  free_array_float (fft_in);

  return out;
}

vector<float>
ifft (const vector<complex<float>>& in)
{
  vector<float> out ((in.size() - 1) * 2);

  /* ensure memory is SSE-aligned (or other vectorized stuff) */
  float *ifft_in = new_array_float (out.size());
  float *ifft_out = new_array_float (out.size());

  /* complex<float> vector and fft_out have the same layout in memory */
  std::copy (in.begin(), in.end(), reinterpret_cast<complex<float> *> (ifft_in));
  fftsr_float (out.size(), ifft_in, ifft_out);

  std::copy (ifft_out, ifft_out + out.size(), &out[0]);

  free_array_float (ifft_out);
  free_array_float (ifft_in);

  return out;
}


