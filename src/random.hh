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

#ifndef AUDIOWMARK_RANDOM_HH
#define AUDIOWMARK_RANDOM_HH

#include <gcrypt.h>
#include <stdint.h>

#include <vector>
#include <string>
#include <random>

class Random
{
public:
  enum class Stream {
    data_up_down = 1,
    sync_up_down = 2,
    speed_clip = 3,
    mix = 4,
    bit_order = 5,
    frame_position = 6
  };
private:
  gcry_cipher_hd_t           aes_ctr_cipher = nullptr;
  gcry_cipher_hd_t           seed_cipher = nullptr;
  std::vector<uint64_t>      buffer;
  size_t                     buffer_pos = 0;

  std::uniform_real_distribution<double> double_dist;

  void die_on_error (const char *func, gcry_error_t error);
public:
  Random (uint64_t seed, Stream stream);
  ~Random();

  typedef uint64_t result_type;

  result_type
  operator()()
  {
    if (buffer_pos == buffer.size())
      refill_buffer();

    return buffer[buffer_pos++];
  }
  static constexpr result_type
  min()
  {
    return 0;
  }
  static constexpr result_type
  max()
  {
    return UINT64_MAX;
  }
  double
  random_double() /* range [0,1) */
  {
    return double_dist (*this);
  }
  void refill_buffer();
  void seed (uint64_t seed, Stream stream);

  template<class T> void
  shuffle (std::vector<T>& result)
  {
    // Fisher–Yates shuffle
    for (size_t i = 0; i < result.size(); i++)
      {
        const uint64_t random_number = (*this)();

        size_t j = i + random_number % (result.size() - i);
        std::swap (result[i], result[j]);
      }
  }

  static void        set_global_test_key (uint64_t seed);
  static void        load_global_key (const std::string& key_file);
  static std::string gen_key();
  static uint64_t    seed_from_hash (const std::vector<float>& floats);
};

#endif /* AUDIOWMARK_RANDOM_HH */
