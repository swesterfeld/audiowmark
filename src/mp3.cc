#include "mp3.hh"

#include <mpg123.h>
#include <assert.h>
#include <stdio.h>
#include <vector>

using std::vector;
using std::string;

bool
mp3_try_load (const string& filename, WavData& wav_data)
{
  int err = mpg123_init();
  assert (err == MPG123_OK);

  mpg123_handle *mh = mpg123_new (NULL, &err);
  assert (err == MPG123_OK);

  mpg123_param (mh, MPG123_ADD_FLAGS, MPG123_QUIET, 0);
  long rate;
  int channels;
  int encoding;
  //
  const long *rates;
  size_t rate_count;
  mpg123_format_none (mh);
  mpg123_rates(&rates, &rate_count);
  for(size_t i=0; i<rate_count; ++i)
    {
      err = mpg123_format (mh, rates[i], MPG123_MONO|MPG123_STEREO, MPG123_ENC_FLOAT_32);
      assert (err == 0);
    }
  printf ("# format: %d\n", err);

  err = mpg123_open (mh, filename.c_str());
  assert (err == 0);
  err = mpg123_getformat (mh, &rate, &channels, &encoding);
  assert (err == 0);
  printf ("# %d\n", err);
  printf ("# %ld %d %d\n", rate, channels, encoding);
  assert (err == MPG123_OK);

  /* ensure that the format will not change */
  mpg123_format_none (mh);
  mpg123_format (mh, rate, channels, encoding);

  printf ("# %zd\n", mpg123_outblock (mh));
  vector<unsigned char> buffer (mpg123_outblock (mh));

  size_t done = 0;
  do
    {
      err = mpg123_read( mh, &buffer[0], buffer.size(), &done );
      assert (err == 0 || err == MPG123_DONE);
      printf ("# done=%zd err=%d\n", done, err);

      float *f = reinterpret_cast<float *> (&buffer[0]);
      for (int i = 0; i < buffer.size() / 8; i++)
        {
          printf ("%f\n", f[i*2]);
        }
    }
  while (done);

  return true;
}
