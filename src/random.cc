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

#include "random.hh"
#include "utils.hh"

#include <regex>

#include <assert.h>

using std::string;
using std::vector;
using std::regex;
using std::regex_match;

static void
gcrypt_init()
{
  static bool init_ok = false;

  if (!init_ok)
    {
      /* version check: start libgcrypt initialization */
      if (!gcry_check_version (GCRYPT_VERSION))
        {
          error ("audiowmark: libgcrypt version mismatch\n");
          exit (1);
        }

      /* disable secure memory (assume we run in a controlled environment) */
      gcry_control (GCRYCTL_DISABLE_SECMEM, 0);

      /* tell libgcrypt that initialization has completed */
      gcry_control (GCRYCTL_INITIALIZATION_FINISHED, 0);

      init_ok = true;
    }
}


static vector<unsigned char> aes_key (16); // 128 bits
static constexpr auto        GCRY_CIPHER = GCRY_CIPHER_AES128;

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

#if 0 /* debugging only */
static void
print (const string& label, const vector<unsigned char>& data)
{
  printf ("%s: ", label.c_str());
  for (auto ch : data)
    printf ("%02x ", ch);
  printf ("\n");
}
#endif

Random::Random (uint64_t start_seed, Stream stream)
{
  gcrypt_init();

  gcry_error_t gcry_ret = gcry_cipher_open (&aes_ctr_cipher, GCRY_CIPHER, GCRY_CIPHER_MODE_CTR, 0);
  die_on_error ("gcry_cipher_open", gcry_ret);

  gcry_ret = gcry_cipher_setkey (aes_ctr_cipher, &aes_key[0], aes_key.size());
  die_on_error ("gcry_cipher_setkey", gcry_ret);

  gcry_ret = gcry_cipher_open (&seed_cipher, GCRY_CIPHER, GCRY_CIPHER_MODE_ECB, 0);
  die_on_error ("gcry_cipher_open", gcry_ret);

  gcry_ret = gcry_cipher_setkey (seed_cipher, &aes_key[0], aes_key.size());
  die_on_error ("gcry_cipher_setkey", gcry_ret);

  seed (start_seed, stream);
}

void
Random::seed (uint64_t seed, Stream stream)
{
  buffer_pos = 0;
  buffer.clear();

  unsigned char plain_text[aes_key.size()];
  unsigned char cipher_text[aes_key.size()];

  memset (plain_text, 0, sizeof (plain_text));
  uint64_to_buffer (seed, &plain_text[0]);

  plain_text[8] = uint8_t (stream);

  gcry_error_t gcry_ret = gcry_cipher_encrypt (seed_cipher, &cipher_text[0], aes_key.size(),
                                                            &plain_text[0],  aes_key.size());
  die_on_error ("gcry_cipher_encrypt", gcry_ret);

  gcry_ret = gcry_cipher_setctr (aes_ctr_cipher, &cipher_text[0], aes_key.size());
  die_on_error ("gcry_cipher_setctr", gcry_ret);
}

Random::~Random()
{
  gcry_cipher_close (aes_ctr_cipher);
  gcry_cipher_close (seed_cipher);
}

void
Random::refill_buffer()
{
  const size_t block_size = 256;
  static unsigned char zeros[block_size] = { 0, };
  unsigned char cipher_text[block_size];

  gcry_error_t gcry_ret = gcry_cipher_encrypt (aes_ctr_cipher, cipher_text, block_size, zeros, block_size);
  die_on_error ("gcry_cipher_encrypt", gcry_ret);

  // print ("AES OUT", {cipher_text, cipher_text + block_size});

  buffer.clear();
  for (size_t i = 0; i < block_size; i += 8)
    buffer.push_back (uint64_from_buffer (cipher_text + i));

  buffer_pos = 0;
}

void
Random::die_on_error (const char *func, gcry_error_t err)
{
  if (err)
    {
      error ("%s failed: %s/%s\n", func, gcry_strsource (err), gcry_strerror (err));

      exit (1); /* can't recover here */
    }
}

void
Random::set_global_test_key (uint64_t key)
{
  uint64_to_buffer (key, &aes_key[0]);
}

void
Random::load_global_key (const string& key_file)
{
  FILE *f = fopen (key_file.c_str(), "r");
  if (!f)
    {
      error ("audiowmark: error opening key file: '%s'\n", key_file.c_str());
      exit (1);
    }

  const regex blank_re (R"(\s*(#.*)?[\r\n]+)");
  const regex key_re (R"(\s*key\s+([0-9a-f]+)\s*(#.*)?[\r\n]+)");

  char buffer[1024];
  int line = 1;
  int keys = 0;
  while (fgets (buffer, 1024, f))
    {
      string s = buffer;

      std::smatch match;
      if (regex_match (s, blank_re))
        {
          /* blank line or comment */
        }
      else if (regex_match (s, match, key_re))
        {
          /* line containing aes key */
          vector<unsigned char> key = hex_str_to_vec (match[1].str());
          if (key.size() != aes_key.size())
            {
              error ("audiowmark: wrong key length in key file '%s', line %d\n => required key length is %zd bits\n", key_file.c_str(), line, aes_key.size() * 8);
              exit (1);
            }
          aes_key = key;
          keys++;
        }
      else
        {
          error ("audiowmark: parse error in key file '%s', line %d\n", key_file.c_str(), line);
          exit (1);
        }
      line++;
    }
  fclose (f);

  if (keys > 1)
    {
      error ("audiowmark: key file '%s' contains more than one key\n", key_file.c_str());
      exit (1);
    }
  if (keys == 0)
    {
      error ("audiowmark: key file '%s' contains no key\n", key_file.c_str());
      exit (1);
    }
}

string
Random::gen_key()
{
  gcrypt_init();

  vector<unsigned char> key (16);
  gcry_randomize (&key[0], 16, /* long term key material strength */ GCRY_VERY_STRONG_RANDOM);
  return vec_to_hex_str (key);
}

uint64_t
Random::seed_from_hash (const vector<float>& floats)
{
  unsigned char hash[20];
  gcry_md_hash_buffer (GCRY_MD_SHA1, hash, &floats[0], floats.size() * sizeof (float));
  return uint64_from_buffer (hash);
}
