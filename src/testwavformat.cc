/*
 * Copyright (C) 2018-2024 Stefan Westerfeld
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

#include <string>
#include <vector>
#include <map>
#include <sndfile.h>
#include <assert.h>
#include <math.h>
#include <cstring>

#include "sfinputstream.hh"
#include "stdoutwavoutputstream.hh"
#include "utils.hh"

using std::string;
using std::vector;

int
main (int argc, char **argv)
{
  std::map<std::string, int> formats
    {
      { "pcm_8",  SF_FORMAT_PCM_U8 },
      { "pcm_16", SF_FORMAT_PCM_16 },
      { "pcm_24", SF_FORMAT_PCM_24 },
      { "pcm_32", SF_FORMAT_PCM_32 },
      { "float",  SF_FORMAT_FLOAT },
      { "double", SF_FORMAT_DOUBLE }
    };

  if (argc == 2 && !strcmp (argv[1], "list"))
    {
      for (auto fi : formats)
        printf ("%s\n", fi.first.c_str());
    }
  else if (argc == 3 && !strcmp (argv[1], "detect"))
    {
      SF_INFO sfinfo = { 0, };
      auto sndfile = sf_open (argv[2], SFM_READ, &sfinfo);
      assert (sndfile);

      for (auto fi : formats)
        {
          if ((fi.second | SF_FORMAT_WAV) == sfinfo.format)
            {
              printf ("%s\n", fi.first.c_str());
              return 0;
            }
        }
      fprintf (stderr, "unsupported format %d\n", sfinfo.format);
      return 1;
    }
  else if (argc == 5 && !strcmp (argv[1], "convert"))
    {
      SFInputStream in;

      std::string in_filename = argv[2];
      std::string out_filename = argv[3];
      std::string out_format = argv[4];

      Error err = in.open (in_filename.c_str());
      if (err)
        {
          fprintf (stderr, "testwavformat: open input failed: %s\n", err.message());
          return 1;
        }
      SF_INFO sfinfo = {0,};

      sfinfo.samplerate = in.sample_rate();
      sfinfo.channels   = in.n_channels();


      sfinfo.format = formats[out_format];
      if (!sfinfo.format)
        {
          fprintf (stderr, "testwavformat: unsupported output format %s\n", out_format.c_str());
          return 1;
        }
      sfinfo.format |= SF_FORMAT_WAV;

      auto sndfile = sf_open (out_filename.c_str(), SFM_WRITE, &sfinfo);
      int error = sf_error (sndfile);
      if (error)
        {
          fprintf (stderr, "%s\n", sf_strerror (sndfile));
          return 1;
        }
      vector<float> samples;
      do
        {
          in.read_frames (samples, 1024);
          sf_count_t count = sf_write_float (sndfile, samples.data(), samples.size());
          assert ((uint64_t) count == samples.size());
        }
      while (samples.size());

      sf_close (sndfile);
    }
  else
    {
      fprintf (stderr, "usage: testwavformat convert <in_filename> <out_filename> <format>\n");
      fprintf (stderr, "or     testwavformat detect <in_filename>\n");
      fprintf (stderr, "or     testwavformat list\n");
      return 1;
    }
}
