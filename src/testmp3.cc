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

#include "mp3inputstream.hh"
#include "wavdata.hh"

using std::string;

int
main (int argc, char **argv)
{
  WavData wd;
  if (argc >= 2)
    {
      if (MP3InputStream::detect (argv[1]))
        {
          MP3InputStream m3i;
          Error err = m3i.open (argv[1]);
          if (err)
            {
              printf ("mp3 open %s failed: %s\n", argv[1], err.message());
              return 1;
            }

          err = wd.load (&m3i);
          if (!err)
            {
              int sec = wd.n_values() / wd.n_channels() / wd.sample_rate();

              printf ("loaded mp3 %s: %d:%02d\n", argv[1], sec / 60, sec % 60);
              if (argc == 3)
                {
                  wd.save (argv[2]);
                  printf ("saved wav: %s\n", argv[2]);
                }
            }
          else
            {
              printf ("mp3 load %s failed: %s\n", argv[1], err.message());
              return 1;
            }
        }
      else
        {
          printf ("mp3 detect %s failed\n", argv[1]);
          return 1;
        }
    }
}
