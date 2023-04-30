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
#include <map>
#include <sstream>

#include <sys/time.h>
#include <assert.h>
#include <stdint.h>

#include "shortcode.hh"

using std::vector;
using std::string;
using std::map;

static double
gettime()
{
  timeval tv;
  gettimeofday (&tv, 0);

  return tv.tv_sec + tv.tv_usec / 1000000.0;
}

vector<int>
generate_error_vector (size_t n, int errors)
{
  vector<int> ev (n);

  while (errors)
    {
      size_t pos = rand() % ev.size();
      if (ev[pos] != 1)
        {
          ev[pos] = 1;
          errors--;
        }
    }
  return ev;
}

int
hamming_weight (const vector<int>& bits)
{
  int w = 0;
  for (auto b : bits)
    w += b;
  return w;
}

double
factorial (int x)
{
  double p = 1;
  for (int i = 1; i <= x; i++)
    p *= i;
  return p;
}

string
number_format (double d)
{
  std::ostringstream buff;
  buff.imbue (std::locale(""));
  buff << (uint64_t) d;
  return buff.str();
}

int
main (int argc, char **argv)
{
  srand (time (NULL));

  if (argc < 2)
    {
      printf ("first argument must be code size (12, 16, 20)\n");
      return 1;
    }
  size_t K = atoi (argv[1]);
  size_t N = short_code_init (K);
  if (!N)
    {
      printf ("bad code size\n");
      return 1;
    }
  printf ("using (%zd,%zd) code\n", N, K);

  if (argc == 2)
    {
      vector<int> in_bits;
      while (in_bits.size() != K)
        in_bits.push_back (rand() & 1);

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
  if (argc == 3 && string (argv[2]) == "perf")
    {
      const double start_t = gettime();
      const size_t runs = 100;
      for (size_t i = 0; i < runs; i++)
        {
          vector<int> in_bits;
          while (in_bits.size() != K)
            in_bits.push_back (rand() & 1);

          vector<int> out_bits = short_decode_blk (short_encode_blk (in_bits));
          assert (out_bits == in_bits);
        }
      printf ("%.1f ms/block\n", (gettime() - start_t) / runs * 1000.0);
    }
  if (argc == 3 && string (argv[2]) == "table")
    {
      map<vector<int>, vector<int>> table;
      vector<int> weight (N + 1);
      for (size_t i = 0; i < size_t (1 << K); i++)
        {
          vector<int> in;
          for (size_t bit = 0; bit < K; bit++)
            {
              if (i & (1 << bit))
                in.push_back (1);
              else
                in.push_back (0);
            }
          vector<int> coded_bits = short_encode_blk (in);
          table[coded_bits] = in;
          weight[hamming_weight (coded_bits)]++;
          printf ("T: ");
          for (auto b : coded_bits)
            printf ("%d", b);
          printf ("\n");

        }
      for (size_t i = 0; i <= N; i++)
        {
          if (weight[i])
            printf ("W %3zd %6d %20s\n", i, weight[i], number_format (factorial (N) / (factorial (i) * factorial (N - i)) / weight[i]).c_str());
        }
      /* decoding test */
      for (auto it : table)
        {
          assert (short_decode_blk (it.first) == it.second);
        }

      const size_t runs = 50LL * 1000 * 1000 * 1000;
      size_t match = 0;
      for (size_t i = 0; i < runs; i++)
        {
          vector<int> in_bits = generate_error_vector (N, K);
          auto it = table.find (in_bits);

          if (it != table.end())
            match++;

          if ((i % 1000000) == 0)
            {
              printf ("%zd / %zd\r", match, i);
              fflush (stdout);
            }
        }
    }
  if (argc == 3 && string (argv[2]) == "distance")
    {
      vector<vector<int>> cwords;
      for (size_t i = 0; i < size_t (1 << K); i++)
        {
          vector<int> in;
          for (size_t bit = 0; bit < K; bit++)
            {
              if (i & (1 << bit))
                in.push_back (1);
              else
                in.push_back (0);
            }
          cwords.push_back (short_encode_blk (in));
        }
      int mhd = 100000;
      for (size_t a = 0; a < cwords.size(); a++)
        {
          for (size_t b = 0; b < cwords.size(); b++)
            {
              if (a != b)
                {
                  int hd = 0;

                  for (size_t c = 0; c < cwords[a].size(); c++)
                    hd += cwords[a][c] ^ cwords[b][c];
                  if (hd < mhd)
                    mhd = hd;
                }
            }
          if ((a & 255) == 0)
            {
              printf ("%zd\r", a);
              fflush (stdout);
            }
        }
      printf ("\n");
      printf ("%d\n", mhd);
    }
}
