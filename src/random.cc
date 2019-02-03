#include "random.hh"
#include "utils.hh"

#include <regex>

#include <assert.h>

using std::string;
using std::vector;
using std::regex;
using std::regex_match;

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

Random::Random (uint64_t seed, Stream stream)
{
  vector<unsigned char> ctr = get_start_counter (seed, stream);

  // print ("CTR", ctr);

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

vector<unsigned char>
Random::get_start_counter (uint64_t seed, Stream stream)
{
  gcry_error_t     gcry_ret;
  gcry_cipher_hd_t cipher_hd;

  gcry_ret = gcry_cipher_open (&cipher_hd, GCRY_CIPHER, GCRY_CIPHER_MODE_ECB, 0);
  die_on_error ("gcry_cipher_open", gcry_ret);

  gcry_ret = gcry_cipher_setkey (cipher_hd, &aes_key[0], aes_key.size());
  die_on_error ("gcry_cipher_setkey", gcry_ret);

  vector<unsigned char> cipher_text (16);
  vector<unsigned char> plain_text (16);

  uint64_to_buffer (seed, &plain_text[0]);

  plain_text[8] = uint8_t (stream);

  // print ("SEED", plain_text);

  gcry_ret = gcry_cipher_encrypt (cipher_hd, &cipher_text[0], cipher_text.size(),
                                             &plain_text[0],  plain_text.size());
  die_on_error ("gcry_cipher_encrypt", gcry_ret);

  gcry_cipher_close (cipher_hd);

  return cipher_text;
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
Random::die_on_error (const char *func, gcry_error_t error)
{
  if (error)
    {
      fprintf (stderr, "%s failed: %s/%s\n", func, gcry_strsource (error), gcry_strerror (error));

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
      fprintf (stderr, "audiowmark: error opening key file: '%s'\n", key_file.c_str());
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
              fprintf (stderr, "audiowmark: wrong key length in key file '%s', line %d\n => required key length is %zd bits\n", key_file.c_str(), line, aes_key.size() * 8);
              exit (1);
            }
          aes_key = key;
          keys++;
        }
      else
        {
          fprintf (stderr, "audiowmark: parse error in key file '%s', line %d\n", key_file.c_str(), line);
          exit (1);
        }
      line++;
    }
  fclose (f);

  if (keys > 1)
    {
      fprintf (stderr, "audiowmark: key file '%s' contains more than one key\n", key_file.c_str());
      exit (1);
    }
  if (keys == 0)
    {
      fprintf (stderr, "audiowmark: key file '%s' contains no key\n", key_file.c_str());
      exit (1);
    }
}

string
Random::gen_key()
{
  vector<unsigned char> key (16);
  gcry_randomize (&key[0], 16, /* long term key material strength */ GCRY_VERY_STRONG_RANDOM);
  return vec_to_hex_str (key);
}
