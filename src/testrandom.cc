#include "utils.hh"

#include <sys/time.h>
#include <gcrypt.h>

using std::vector;
using std::string;

static double
gettime()
{
  timeval tv;
  gettimeofday (&tv, 0);

  return tv.tv_sec + tv.tv_usec / 1000000.0;
}

static constexpr auto GCRY_CIPHER = GCRY_CIPHER_AES128;

class Random
{
  vector<unsigned char> aes_key = vector<unsigned char> (16);
  gcry_cipher_hd_t      aes_ctr_cipher;

  void
  die_on_error (const char *func, gcry_error_t error)
  {
    if (error)
      {
        fprintf (stderr, "%s failed: %s/%s\n", func, gcry_strsource (error), gcry_strerror (error));

        exit (1); /* can't recover here */
      }
  }
public:
  enum class Stream {
    up_down = 1,
    mix = 2,
    bit_order = 3
  };
  Random (uint64_t seed, Stream stream)
  {
    vector<unsigned char> ctr = get_start_counter (seed, stream);
    printf ("CTR: ");
    for (auto ch : ctr)
      printf ("%02x ", ch);
    printf ("\n");

    gcry_error_t gcry_ret = gcry_cipher_open (&aes_ctr_cipher, GCRY_CIPHER, GCRY_CIPHER_MODE_CTR, 0);
    die_on_error ("gcry_cipher_open", gcry_ret);

    gcry_ret = gcry_cipher_setkey (aes_ctr_cipher, &aes_key[0], aes_key.size());
    die_on_error ("gcry_cipher_setkey", gcry_ret);

    gcry_ret = gcry_cipher_setctr (aes_ctr_cipher, &ctr[0], ctr.size());
    die_on_error ("gcry_cipher_setctr", gcry_ret);
  }
  ~Random()
  {
    gcry_cipher_close (aes_ctr_cipher);
  }
  vector<unsigned char>
  get_start_counter (uint64_t seed, Stream stream)
  {
    gcry_error_t     gcry_ret;
    gcry_cipher_hd_t cipher_hd;

    gcry_ret = gcry_cipher_open (&cipher_hd, GCRY_CIPHER, GCRY_CIPHER_MODE_ECB, 0);
    die_on_error ("gcry_cipher_open", gcry_ret);

    gcry_ret = gcry_cipher_setkey (cipher_hd, &aes_key[0], aes_key.size());
    die_on_error ("gcry_cipher_setkey", gcry_ret);

    vector<unsigned char> cipher_text (16);
    vector<unsigned char> plain_text (16);

    /* this has to be endian independent: use big endian order */
    plain_text[0] = seed >> 56;
    plain_text[1] = seed >> 48;
    plain_text[2] = seed >> 40;
    plain_text[3] = seed >> 32;
    plain_text[4] = seed >> 24;
    plain_text[5] = seed >> 16;
    plain_text[6] = seed >> 8;
    plain_text[7] = seed;

    plain_text[8] = uint8_t (stream);

    printf ("[[ ");
    for (auto ch : plain_text)
      printf ("%02x ", ch);
    printf (" ]]\n");

    gcry_ret = gcry_cipher_encrypt (cipher_hd, &cipher_text[0], cipher_text.size(),
                                               &plain_text[0],  plain_text.size());
    die_on_error ("gcry_cipher_encrypt", gcry_ret);

    gcry_cipher_close (cipher_hd);

    printf ("[[ ");
    for (auto ch : cipher_text)
      printf ("%02x ", ch);
    printf (" ]]\n");

    return cipher_text;
  }
  uint64_t operator()()
  {
    const size_t block_size = 8;
    unsigned char zeros[block_size] = { 0, };
    unsigned char cipher_text[block_size];

    gcry_error_t gcry_ret = gcry_cipher_encrypt (aes_ctr_cipher, cipher_text, block_size, zeros, block_size);
    die_on_error ("gcry_cipher_encrypt", gcry_ret);

    printf ("[[ ");
    for (auto ch : cipher_text)
      printf ("%02x ", ch);
    printf (" ]]\n");
    /* this has to be endian independent: use big endian order */
    uint64_t result = (uint64_t (cipher_text[0]) << 56)
                    + (uint64_t (cipher_text[1]) << 48)
                    + (uint64_t (cipher_text[2]) << 40)
                    + (uint64_t (cipher_text[3]) << 32)
                    + (uint64_t (cipher_text[4]) << 24)
                    + (uint64_t (cipher_text[5]) << 16)
                    + (uint64_t (cipher_text[6]) << 8)
                    + cipher_text[7];
    return result;
  }
};

int
main (int argc, char **argv)
{
  Random rng (0xf00f1234b00b5678U, Random::Stream::bit_order);
  for (size_t i = 0; i < 20; i++)
    {
      uint64_t x = rng();
      printf ("%016lx\n", x);
    }
}
