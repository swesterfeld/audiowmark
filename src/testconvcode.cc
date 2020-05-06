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
#include "convcode.hh"

#include <random>

#include <assert.h>

using std::vector;
using std::string;

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

static bool
no_case_equal (const string& s1, const string& s2)
{
  if (s1.size() != s2.size())
    return false;

  return std::equal (s1.begin(), s1.end(), s2.begin(),
                     [] (char c1, char c2) -> bool { return tolower (c1) == tolower (c2);});
}

int
main (int argc, char **argv)
{
  string btype = (argc > 1) ? argv[1] : "";
  ConvBlockType block_type;

  if (no_case_equal (btype, "A"))
    block_type = ConvBlockType::a;
  else if (no_case_equal (btype, "B"))
    block_type = ConvBlockType::b;
  else if (no_case_equal (btype, "AB"))
    block_type = ConvBlockType::ab;
  else
    {
      printf ("first argument must be A, B, or AB\n");
      return 1;
    }

  if (argc == 2)
    {
      vector<int> in_bits = bit_str_to_vec ("80f12381");

      printf ("input vector (k=%zd):  ", in_bits.size());
      for (auto b : in_bits)
        printf ("%d", b);
      printf ("\n");

      vector<int> coded_bits = conv_encode (block_type, in_bits);
      printf ("coded vector (n=%zd): ", coded_bits.size());
      for (auto b : coded_bits)
        printf ("%d", b);
      printf ("\n");
      printf ("coded hex: %s\n", bit_vec_to_str (coded_bits).c_str());

      assert (coded_bits.size() == conv_code_size (block_type, in_bits.size()));

      vector<int> decoded_bits = conv_decode_hard (block_type, coded_bits);
      printf ("output vector (k=%zd): ", decoded_bits.size());
      for (auto b : decoded_bits)
        printf ("%d", b);
      printf ("\n");

      assert (decoded_bits.size() == in_bits.size());
      int errors = 0;
      for (size_t i = 0; i < decoded_bits.size(); i++)
        if (decoded_bits[i] != in_bits[i])
          errors++;
      printf ("decoding errors: %d\n", errors);
    }
  if (argc == 3 && string (argv[2]) == "error")
    {
      size_t max_bit_errors = conv_code_size (block_type, 128) * 0.5;

      for (size_t bit_errors = 0; bit_errors < max_bit_errors; bit_errors++)
        {
          size_t coded_bit_count = 0;
          int bad_decode = 0;
          constexpr int test_size = 20;

          for (int i = 0; i < test_size; i++)
            {
              vector<int> in_bits;
              while (in_bits.size() != 128)
                in_bits.push_back (rand() & 1);

              vector<int> coded_bits = conv_encode (block_type, in_bits);
              coded_bit_count = coded_bits.size();

              vector<int> error_bits = generate_error_vector (coded_bits.size(), bit_errors);
              for (size_t pos = 0; pos < coded_bits.size(); pos++)
                coded_bits[pos] ^= error_bits[pos];

              vector<int> decoded_bits = conv_decode_hard (block_type, coded_bits);

              assert (decoded_bits.size() == 128);

              int errors = 0;
              for (size_t i = 0; i < 128; i++)
                if (decoded_bits[i] != in_bits[i])
                  errors++;
              if (errors > 0)
                bad_decode++;
            }
          printf ("%f %f\n", (100.0 * bit_errors) / coded_bit_count, (100.0 * bad_decode) / test_size);
        }
    }
  if (argc == 3 && string (argv[2]) == "soft-error")
    {
      for (double stddev = 0; stddev < 1.5; stddev += 0.01)
        {
          size_t coded_bit_count = 0;
          int bad_decode1 = 0, bad_decode2 = 0;
          constexpr int test_size = 20;

          int local_be = 0;
          for (int i = 0; i < test_size; i++)
            {
              vector<int> in_bits;
              while (in_bits.size() != 128)
                in_bits.push_back (rand() & 1);

              vector<int> coded_bits = conv_encode (block_type, in_bits);
              coded_bit_count = coded_bits.size();

              std::default_random_engine generator;
              std::normal_distribution<double> dist (0, stddev);

              vector<float> recv_bits;
              for (auto b : coded_bits)
                recv_bits.push_back (b + dist (generator));

              vector<int> decoded_bits1 = conv_decode_soft (block_type, recv_bits);

              vector<int> recv_hard_bits;
              for (auto b : recv_bits)
                recv_hard_bits.push_back ((b > 0.5) ? 1 : 0);

              for (size_t x = 0; x < recv_hard_bits.size(); x++)
                local_be += coded_bits[x] ^ recv_hard_bits[x];

              vector<int> decoded_bits2 = conv_decode_hard (block_type, recv_hard_bits);

              assert (decoded_bits1.size() == 128);
              assert (decoded_bits2.size() == 128);

              int e1 = 0;
              int e2 = 0;
              for (size_t i = 0; i < 128; i++)
                {
                  if (decoded_bits1[i] != in_bits[i])
                    e1++;
                  if (decoded_bits2[i] != in_bits[i])
                    e2++;
                }
              if (e1)
                bad_decode1++;
              if (e2)
                bad_decode2++;
            }
          printf ("%f %f %f\n", double (100 * local_be) / test_size / coded_bit_count, (100.0 * bad_decode1) / test_size, (100.0 * bad_decode2) / test_size);
        }
    }
  if (argc == 3 && string (argv[2]) == "perf")
    {
      vector<int> in_bits;
      while (in_bits.size() != 128)
        in_bits.push_back (rand() & 1);

      const double start_t = get_time();
      const size_t runs = 20;
      for (size_t i = 0; i < runs; i++)
        {
          vector<int> out_bits = conv_decode_hard (block_type, conv_encode (block_type, in_bits));
          assert (out_bits == in_bits);
        }
      printf ("%.1f ms/block\n", (get_time() - start_t) / runs * 1000.0);
    }
  if (argc == 3 && string (argv[2]) == "table")
    conv_print_table (block_type);
}
