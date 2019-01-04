#include "random.hh"

using std::string;

static std::vector<unsigned char> aes_key (16);
static constexpr auto             GCRY_CIPHER = GCRY_CIPHER_AES128;

static void
uint64_to_buffer (uint64_t       u,
                  unsigned char *buffer)
{
  /* this has to be endian independent: use big endian order */
  buffer[0] = u >> 56;
  buffer[1] = u >> 48;
  buffer[2] = u >> 40;
  buffer[3] = u >> 32;
  buffer[4] = u >> 24;
  buffer[5] = u >> 16;
  buffer[6] = u >> 8;
  buffer[7] = u;
}

static uint64_t
uint64_from_buffer (unsigned char *buffer)
{
  /* this has to be endian independent: use big endian order */
  return (uint64_t (buffer[0]) << 56)
       + (uint64_t (buffer[1]) << 48)
       + (uint64_t (buffer[2]) << 40)
       + (uint64_t (buffer[3]) << 32)
       + (uint64_t (buffer[4]) << 24)
       + (uint64_t (buffer[5]) << 16)
       + (uint64_t (buffer[6]) << 8)
       + buffer[7];
}

Random::Random (uint64_t seed, Stream stream)
{
  std::vector<unsigned char> ctr = get_start_counter (seed, stream);
#if 0
  printf ("CTR: ");
  for (auto ch : ctr)
    printf ("%02x ", ch);
  printf ("\n");
#endif

  gcry_error_t gcry_ret = gcry_cipher_open (&aes_ctr_cipher, GCRY_CIPHER, GCRY_CIPHER_MODE_CTR, 0);
  die_on_error ("gcry_cipher_open", gcry_ret);

  gcry_ret = gcry_cipher_setkey (aes_ctr_cipher, &aes_key[0], aes_key.size());
  die_on_error ("gcry_cipher_setkey", gcry_ret);

  gcry_ret = gcry_cipher_setctr (aes_ctr_cipher, &ctr[0], ctr.size());
  die_on_error ("gcry_cipher_setctr", gcry_ret);
}

Random::~Random()
{
  gcry_cipher_close (aes_ctr_cipher);
}

std::vector<unsigned char>
Random::get_start_counter (uint64_t seed, Stream stream)
{
  gcry_error_t     gcry_ret;
  gcry_cipher_hd_t cipher_hd;

  gcry_ret = gcry_cipher_open (&cipher_hd, GCRY_CIPHER, GCRY_CIPHER_MODE_ECB, 0);
  die_on_error ("gcry_cipher_open", gcry_ret);

  gcry_ret = gcry_cipher_setkey (cipher_hd, &aes_key[0], aes_key.size());
  die_on_error ("gcry_cipher_setkey", gcry_ret);

  std::vector<unsigned char> cipher_text (16);
  std::vector<unsigned char> plain_text (16);

  uint64_to_buffer (seed, &plain_text[0]);

  plain_text[8] = uint8_t (stream);

#if 0
  printf ("[[ ");
  for (auto ch : plain_text)
    printf ("%02x ", ch);
  printf (" ]]\n");
#endif

  gcry_ret = gcry_cipher_encrypt (cipher_hd, &cipher_text[0], cipher_text.size(),
                                             &plain_text[0],  plain_text.size());
  die_on_error ("gcry_cipher_encrypt", gcry_ret);

  gcry_cipher_close (cipher_hd);

#if 0
  printf ("[[ ");
  for (auto ch : cipher_text)
    printf ("%02x ", ch);
  printf (" ]]\n");
#endif

  return cipher_text;
}

uint64_t
Random::operator()()
{
  const size_t block_size = 8;
  unsigned char zeros[block_size] = { 0, };
  unsigned char cipher_text[block_size];

  gcry_error_t gcry_ret = gcry_cipher_encrypt (aes_ctr_cipher, cipher_text, block_size, zeros, block_size);
  die_on_error ("gcry_cipher_encrypt", gcry_ret);

#if 0
  printf ("[[ ");
  for (auto ch : cipher_text)
    printf ("%02x ", ch);
  printf (" ]]\n");
#endif
  return uint64_from_buffer (&cipher_text[0]);
}

void
Random::set_global_test_key (uint64_t key)
{
  uint64_to_buffer (key, &aes_key[0]);
}

std::string
Random::gen_key()
{
  unsigned char key[16];
  gcry_randomize (key, 16, /* long term key material strength */ GCRY_VERY_STRONG_RANDOM);
  string s;
  for (auto k : key)
    {
      char buffer[256];

      sprintf (buffer, "%02x", k);
      s += buffer;
    }
  return s;
}
