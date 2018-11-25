#include <string.h>
#include <math.h>
#include <string>

#include <fftw3.h>

#include "wavdata.hh"

using std::string;
using std::vector;

namespace Params
{
  static constexpr int frame_size = 1024;
}

inline double
window_cos (double x) /* von Hann window */
{
  if (fabs (x) > 1)
    return 0;
  return 0.5 * cos (x * M_PI) + 0.5;
}

inline double
window_hamming (double x) /* sharp (rectangle) cutoffs at boundaries */
{
  if (fabs (x) > 1)
    return 0;

  return 0.54 + 0.46 * cos (M_PI * x);
}

int
frame_count (WavData& wav_data)
{
  return (wav_data.n_values() / wav_data.n_channels() + (Params::frame_size - 1)) / Params::frame_size;
}

/*
 * get one audio frame, Params::frame_size samples if available
 *
 * in case of stereo: deinterleave
 */
vector<float>
get_frame (WavData& wav_data, int f, int ch)
{
  auto& samples = wav_data.samples();

  vector<float> result;

  size_t pos = (f * Params::frame_size) * wav_data.n_channels() + ch;
  for (size_t x = 0; x < Params::frame_size; x++)
    {
      if (pos < samples.size())
        result.push_back (samples[pos]);

      pos += wav_data.n_channels();
    }
  return result;
}

float *
new_array_float (size_t N)
{
  const size_t N_2 = N + 2; /* extra space for r2c extra complex output */

  return (float *) fftwf_malloc (sizeof (float) * N_2);
}

float *
free_array_float (float *f)
{
  fftwf_free (f);
}

void
fftar_float (size_t N, float *in, float *out)
{
  static fftwf_plan plan = nullptr; // FIXME: should be one plan per fft size

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
  static fftwf_plan plan = nullptr; // FIXME: should be one plan per fft size

  if (!plan)
    {
      float *plan_in = new_array_float (N);
      float *plan_out = new_array_float (N);
      plan = fftwf_plan_dft_c2r_1d (N, (fftwf_complex *) plan_in, plan_out, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);

      // we add code for saving plans here, and use patient planning
    }
  fftwf_execute_dft_c2r (plan, (fftwf_complex *)in, out);
}

int
add_watermark (const string& infile, const string& outfile, const string& bits)
{
  printf ("loading %s\n", infile.c_str());

  WavData wav_data;
  if (!wav_data.load (infile))
    {
      fprintf (stderr, "audiowmark: error loading %s: %s\n", infile.c_str(), wav_data.error_blurb());
      return 1;
    }
  printf ("channels: %d, samples: %zd, mix_freq: %f\n", wav_data.n_channels(), wav_data.n_values(), wav_data.mix_freq());

  for (int f = 0; f < frame_count (wav_data); f++)
    {
      for (int ch = 0; ch < wav_data.n_channels(); ch++)
        {
          vector<float> frame = get_frame (wav_data, f, ch);
          if (frame.size() == Params::frame_size)
            {
              /* windowing */
              double window_weight = 0;
              for (size_t i = 0; i < frame.size(); i++)
                {
                  const double fsize_2 = frame.size() / 2.0;
                  // const double win =  window_cos ((i - fsize_2) / fsize_2);
                  const double win = window_hamming ((i - fsize_2) / fsize_2);
                  //const double win = 1;
                  frame[i] *= win;
                  window_weight += win;
                }

              /* to get normalized fft output corrected by window weight */
              for (size_t i = 0; i < frame.size(); i++)
                frame[i] *= 2.0 / window_weight;

              /* FFT transform */
              float *fft_in = new_array_float (frame.size());
              float *fft_out = new_array_float (frame.size());
              std::copy (frame.begin(), frame.end(), fft_in);
              fftar_float (frame.size(), fft_in, fft_out);

              for (size_t i = 0; i <= Params::frame_size / 2; i++)
                {
                  const double re = fft_out[i * 2];
                  const double im = fft_out[i * 2 + 1];
                  const double mag = sqrt (re * re + im * im);
                  printf ("fft %d %d %zd %f\n", f, ch, i, mag);
                }

              free_array_float (fft_out);
              free_array_float (fft_in);
            }
        }
    }

  if (!wav_data.save (outfile))
    {
      fprintf (stderr, "audiowmark: error saving %s: %s\n", outfile.c_str(), wav_data.error_blurb());
      return 1;
    }
}

int
main (int argc, char **argv)
{
  if (strcmp (argv[1], "add") == 0 && argc == 5)
    {
      return add_watermark (argv[2], argv[3], argv[4]);
    }
  else
    {
      fprintf (stderr, "audiowmark: error parsing commandline args\n");
      return 1;
    }
}
