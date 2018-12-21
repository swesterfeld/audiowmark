#include "utils.hh"
#include "convcode.hh"

#include <random>

#include <assert.h>

using std::vector;
using std::string;

vector<int>
generate_error_vector (size_t n, int errors)
{
  vector<int> ev (n);

  while (errors)
    {
      size_t pos = rand() % ev.size();
      if (ev[pos] != 1)
        {
          ev[pos] = 1;
          errors--;
        }
    }
  return ev;
}
int
main (int argc, char **argv)
{
  if (argc == 1)
    {
      vector<int> in_bits = bit_str_to_vec ("80f12381");

      printf ("input vector (k=%zd):  ", in_bits.size());
      for (auto b : in_bits)
        printf ("%d", b);
      printf ("\n");

      vector<int> coded_bits = conv_encode (in_bits);
      printf ("coded vector (n=%zd): ", coded_bits.size());
      for (auto b : coded_bits)
        printf ("%d", b);
      printf ("\n");
      printf ("coded hex: %s\n", bit_vec_to_str (coded_bits).c_str());

      assert (coded_bits.size() == conv_code_size (in_bits.size()));

      vector<int> decoded_bits = conv_decode_hard (coded_bits);
      printf ("output vector (k=%zd): ", decoded_bits.size());
      for (auto b : decoded_bits)
        printf ("%d", b);
      printf ("\n");

      assert (decoded_bits.size() == in_bits.size());
      int errors = 0;
      for (size_t i = 0; i < decoded_bits.size(); i++)
        if (decoded_bits[i] != in_bits[i])
          errors++;
      printf ("decoding errors: %d\n", errors);
    }
  if (argc == 2 && string (argv[1]) == "error")
    {
      for (size_t bit_errors = 0; bit_errors < 100; bit_errors++)
        {
          size_t coded_bit_count = 0;
          int bad_decode = 0;
          constexpr int test_size = 20;

          for (int i = 0; i < test_size; i++)
            {
              vector<int> in_bits;
              while (in_bits.size() != 128)
                in_bits.push_back (rand() & 1);

              vector<int> coded_bits = conv_encode (in_bits);
              coded_bit_count = coded_bits.size();

              vector<int> error_bits = generate_error_vector (coded_bits.size(), bit_errors);
              for (size_t pos = 0; pos < coded_bits.size(); pos++)
                coded_bits[pos] ^= error_bits[pos];

              vector<int> decoded_bits = conv_decode_hard (coded_bits);

              assert (decoded_bits.size() == 128);

              int errors = 0;
              for (size_t i = 0; i < 128; i++)
                if (decoded_bits[i] != in_bits[i])
                  errors++;
              if (errors > 0)
                bad_decode++;
            }
          printf ("%f %f\n", (100.0 * bit_errors) / coded_bit_count, (100.0 * bad_decode) / test_size);
        }
    }
  if (argc == 2 && string (argv[1]) == "soft-error")
    {
      for (double stddev = 0; stddev < 1; stddev += 0.01)
        {
          size_t coded_bit_count = 0;
          int bad_decode1 = 0, bad_decode2 = 0;
          constexpr int test_size = 20;

          int local_be = 0;
          for (int i = 0; i < test_size; i++)
            {
              vector<int> in_bits;
              while (in_bits.size() != 128)
                in_bits.push_back (rand() & 1);

              vector<int> coded_bits = conv_encode (in_bits);
              coded_bit_count = coded_bits.size();

              std::default_random_engine generator;
              std::normal_distribution<double> dist (0, stddev);

              vector<float> recv_bits;
              for (auto b : coded_bits)
                recv_bits.push_back (b + dist (generator));

              vector<int> decoded_bits1 = conv_decode_soft (recv_bits);

              vector<int> recv_hard_bits;
              for (auto b : recv_bits)
                recv_hard_bits.push_back ((b > 0.5) ? 1 : 0);

              for (size_t x = 0; x < recv_hard_bits.size(); x++)
                local_be += coded_bits[x] ^ recv_hard_bits[x];

              vector<int> decoded_bits2 = conv_decode_hard (recv_hard_bits);

              assert (decoded_bits1.size() == 128);
              assert (decoded_bits2.size() == 128);

              int e1 = 0;
              int e2 = 0;
              for (size_t i = 0; i < 128; i++)
                {
                  if (decoded_bits1[i] != in_bits[i])
                    e1++;
                  if (decoded_bits2[i] != in_bits[i])
                    e2++;
                }
              if (e1)
                bad_decode1++;
              if (e2)
                bad_decode2++;
            }
          printf ("%f %f %f\n", double (100 * local_be) / test_size / coded_bit_count, (100.0 * bad_decode1) / test_size, (100.0 * bad_decode2) / test_size);
        }
    }
}
