#include "random.hh"

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
