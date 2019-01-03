#ifndef AUDIOWMARK_RANDOM_HH
#define AUDIOWMARK_RANDOM_HH

#include <gcrypt.h>
#include <stdint.h>

#include <vector>

class Random
{
public:
  enum class Stream {
    up_down = 1,
    mix = 2,
    bit_order = 3
  };
private:
  std::vector<unsigned char> aes_key = std::vector<unsigned char> (16);
  gcry_cipher_hd_t           aes_ctr_cipher;

  static constexpr auto      GCRY_CIPHER = GCRY_CIPHER_AES128;

  void
  die_on_error (const char *func, gcry_error_t error)
  {
    if (error)
      {
        fprintf (stderr, "%s failed: %s/%s\n", func, gcry_strsource (error), gcry_strerror (error));

        exit (1); /* can't recover here */
      }
  }
  std::vector<unsigned char> get_start_counter (uint64_t seed, Stream stream);
public:
  Random (uint64_t seed, Stream stream);
  ~Random();
  uint64_t operator()();
};

#endif /* AUDIOWMARK_RANDOM_HH */
