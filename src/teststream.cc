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

#include <string>
#include <vector>
#include <sndfile.h>
#include <assert.h>
#include <math.h>

#include "sfinputstream.hh"
#include "stdoutwavoutputstream.hh"
#include "utils.hh"

using std::string;
using std::vector;

int
main (int argc, char **argv)
{
  SFInputStream in;
  StdoutWavOutputStream out;

  std::string filename = (argc >= 2) ? argv[1] : "-";
  Error err = in.open (filename.c_str());
  if (err)
    {
      fprintf (stderr, "teststream: open input failed: %s\n", err.message());
      return 1;
    }
  err = out.open (in.n_channels(), in.sample_rate(), 16, in.n_frames());
  if (err)
    {
      fprintf (stderr, "teststream: open output failed: %s\n", err.message());
      return 1;
    }
  vector<float> samples;
  do
    {
      in.read_frames (samples, 1024);
      out.write_frames (samples);
    }
  while (samples.size());
}
