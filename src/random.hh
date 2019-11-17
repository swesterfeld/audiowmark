#ifndef AUDIOWMARK_RANDOM_HH
#define AUDIOWMARK_RANDOM_HH

#include <gcrypt.h>
#include <stdint.h>

#include <vector>
#include <string>

class Random
{
public:
  enum class Stream {
    data_up_down = 1,
    sync_up_down = 2,
    pad_up_down = 3,
    mix = 4,
    bit_order = 5,
    frame_position = 6
  };
private:
  gcry_cipher_hd_t           aes_ctr_cipher = nullptr;
  gcry_cipher_hd_t           seed_cipher = nullptr;
  std::vector<uint64_t>      buffer;
  size_t                     buffer_pos = 0;

  void die_on_error (const char *func, gcry_error_t error);
public:
  Random (uint64_t seed, Stream stream);
  ~Random();

  uint64_t
  operator()()
  {
    if (buffer_pos == buffer.size())
      refill_buffer();

    return buffer[buffer_pos++];
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
};

#endif /* AUDIOWMARK_RANDOM_HH */
