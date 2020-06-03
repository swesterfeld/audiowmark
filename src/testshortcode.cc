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

#include <vector>

#include <sys/time.h>
#include <assert.h>

#include "shortcode.hh"

using std::vector;
using std::string;

static double
gettime()
{
  timeval tv;
  gettimeofday (&tv, 0);

  return tv.tv_sec + tv.tv_usec / 1000000.0;
}

int
main (int argc, char **argv)
{
  if (argc == 1)
    {
      vector<int> in_bits = { 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0 };

      printf ("in: ");
      for (auto b : in_bits)
        printf ("%d", b);
      printf ("\n");

      printf ("coded: ");
      vector<int> coded_bits = short_encode_blk (in_bits);
      for (auto b : coded_bits)
        printf ("%d", b);
      printf ("\n");

      vector<int> decoded_bits = short_decode_blk (coded_bits);
      printf ("out: ");
      for (auto b : decoded_bits)
        printf ("%d", b);
      printf ("\n");
    }
  if (argc == 2 && string (argv[1]) == "perf")
    {
      vector<int> in_bits;
      while (in_bits.size() != 16)
        in_bits.push_back (rand() & 1);

      const double start_t = gettime();
      const size_t runs = 20;
      for (size_t i = 0; i < runs; i++)
        {
          vector<int> out_bits = short_decode_blk (short_encode_blk (in_bits));
          assert (out_bits == in_bits);
        }
      printf ("%.1f ms/block\n", (gettime() - start_t) / runs * 1000.0);
    }
}
