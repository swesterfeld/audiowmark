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
#include <assert.h>
#include <math.h>
#include <gcrypt.h>

#include <array>
#include <set>

#include "rawconverter.hh"
#include "config.h"

using std::vector;
using std::string;

void
test_int16 (const char *label, const vector<float>& in_samples, Encoding encoding)
{
  RawFormat format;
  format.set_bit_depth (16);
  format.set_encoding (encoding);
#ifdef WORDS_BIGENDDIAN
  format.set_endian (RawFormat::BIG);
#else
  format.set_endian (RawFormat::LITTLE);
#endif

  Error error;
  RawConverter *converter = RawConverter::create (format, error);
  if (error)
    {
      printf ("error: %s\n", error.message());
      exit (1);
    }
  vector<unsigned char> conv (in_samples.size() * 4);
  converter->to_raw (in_samples.data(), conv.data(), in_samples.size());
  float max_diff = 0;
  for (size_t i = 0; i < in_samples.size(); i++)
    {
      if (encoding == Encoding::SIGNED)
        max_diff = std::max (fabs (in_samples[i] * (1 << 15) - ((int16_t *)conv.data())[i]), max_diff);
      else
        max_diff = std::max (fabs ((in_samples[i] + 1) * (1 << 15) - ((uint16_t *)conv.data())[i]), max_diff);
    }
  printf ("%s: max_diff = %f [ should be less than 1.01 ]\n", label, max_diff);
  assert (max_diff < 1.01);
}

string
hash_bytes (vector<unsigned char>& bytes)
{
  string result;
  std::array<unsigned char, 20> hash;
  gcry_md_hash_buffer (GCRY_MD_SHA1, hash.data(), bytes.data(), bytes.size());
  for (auto ch : hash)
    result += string_printf ("%02x", ch);
  return result;
}

int
main (int argc, char **argv)
{
  std::set<string> hashes;
  Error error;
  RawFormat format;
  uint64_t K = 33452759; // prime
  vector<float> in_samples (K), out_samples (K);
  vector<unsigned char> bytes (K * 4);
  for (uint64_t k = 0; k <= K; k++)
    in_samples[k] = (-1 + double (2 * k) / K);

  test_int16 ("int16", in_samples, Encoding::SIGNED);
  test_int16 ("uint16", in_samples, Encoding::UNSIGNED);
  printf ("\n");
  for (auto bit_depth : { 16, 24, 32 })
    {
      for (auto encoding : { Encoding::SIGNED, Encoding::UNSIGNED })
        {
          for (auto endian : { RawFormat::LITTLE, RawFormat::BIG })
            {
              format.set_bit_depth (bit_depth);
              format.set_encoding (encoding);
              format.set_endian (endian);

              RawConverter *converter = RawConverter::create (format, error);
              if (error)
                {
                  printf ("error: %s\n", error.message());
                  return 1;
                }

              std::fill (bytes.begin(), bytes.end(), 0);

              double time1 = get_time();
              converter->to_raw (in_samples.data(), bytes.data(), in_samples.size());
              double time2 = get_time();
              converter->from_raw (bytes.data(), out_samples.data(), in_samples.size());
              double time3 = get_time();

              double max_err = 0;
              for (size_t i = 0; i < in_samples.size(); i++)
                max_err = std::max (max_err, std::abs (double (in_samples[i]) - double (out_samples[i])));
              double ebits = log2 (max_err);
              printf ("%s %d %s endian %f",
                  format.encoding() == Encoding::SIGNED ? "signed" : "unsigned",
                  format.bit_depth(),
                  format.endian() == RawFormat::LITTLE ? "little" : "big",
                  ebits);
              double min_ebits = -format.bit_depth() + 0.9;
              double max_ebits = -format.bit_depth() + 1;
              double ns_per_sample_to = (time2 - time1) * 1e9 / K;
              double ns_per_sample_from = (time3 - time2) * 1e9 / K;
              printf (" (should be in [%.2f,%.2f]) - to raw: %f ns/sample - from raw: %f ns/sample\n",
                  min_ebits, max_ebits,
                  ns_per_sample_to, ns_per_sample_from);
              assert (ebits <= max_ebits && ebits >= min_ebits);

              /* every raw converted buffer should be different */
              string hash = hash_bytes (bytes);
              assert (hashes.count (hash) == 0);
              hashes.insert (hash);
            }
        }
      printf ("\n");
    }
}
