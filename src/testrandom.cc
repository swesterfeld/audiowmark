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

#include "utils.hh"
#include "random.hh"

#include <inttypes.h>

using std::vector;
using std::string;

int
main (int argc, char **argv)
{
  Random rng (0xf00f1234b00b5678U, Random::Stream::bit_order);
  for (size_t i = 0; i < 20; i++)
    {
      uint64_t x = rng();
      printf ("%016" PRIx64 "\n", x);
    }
  for (size_t i = 0; i < 20; i++)
    printf ("%f\n", rng.random_double());

  uint64_t s = 0;
  double t_start = get_time();
  size_t runs = 25000000;
  for (size_t i = 0; i < runs; i++)
    {
      s += rng();
    }
  double t_end = get_time();
  printf ("s=%016" PRIx64 "\n\n", s);

  printf ("%f Mvalues/sec\n", runs / (t_end - t_start) / 1000000);
}
