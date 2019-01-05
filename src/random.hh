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
    up_down = 1,
    mix = 2,
    bit_order = 3
  };
private:
  gcry_cipher_hd_t           aes_ctr_cipher;
  std::vector<uint64_t>      buffer;
  size_t                     buffer_pos = 0;

  std::vector<unsigned char> get_start_counter (uint64_t seed, Stream stream);

  void die_on_error (const char *func, gcry_error_t error);
public:
  Random (uint64_t seed, Stream stream);
  ~Random();
  uint64_t operator()();

  static void        set_global_test_key (uint64_t seed);
  static void        load_global_key (const std::string& key_file);
  static std::string gen_key();
};

#endif /* AUDIOWMARK_RANDOM_HH */
