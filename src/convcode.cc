#include "utils.hh"

#include <array>
#include <assert.h>

using std::vector;

int
parity (unsigned int v)
{
  int p = 0;

  while (v)
    {
      p ^= (v & 1);
      v >>= 1;
    }
  return p;
}

constexpr  unsigned int rate        = 3;
constexpr  unsigned int order       = 9;
constexpr  auto         generators  = std::array<unsigned,3> { 0557, 0663, 0711 };

vector<int>
conv_encode (const vector<int>& in_bits)
{
  vector<int> out_vec;
  vector<int> vec = in_bits;

  /* termination: bring encoder into all-zero state */
  for (unsigned int i = 0; i < order; i++)
    vec.push_back (0);

  unsigned int reg = 0;
  for (auto b : vec)
    {
      reg = (reg << 1) | b;

      for (auto poly : generators)
        {
          int out_bit = parity (reg & poly);

          out_vec.push_back (out_bit);
        }
    }
  return out_vec;
}

vector<int>
conv_decode (const vector<int>& coded_bits)
{
  vector<int> decoded_bits;

  assert (coded_bits.size() % rate == 0);

  int state = 0;

  for (size_t i = 0; i < coded_bits.size(); i += rate)
    {
      int best_delta = order;
      int best_bit = 0;
      for (int bit = 0; bit < 2; bit++)
        {
          int new_state = (state << 1) | bit;
          int delta = 0;

          for (size_t p = 0; p < generators.size(); p++)
            {
              int out_bit = parity (new_state & generators[p]);
              if (out_bit != coded_bits[i + p])
                delta += 1;
            }
          if (delta < best_delta)
            {
              best_bit = bit;
              best_delta = delta;
            }
        }
      decoded_bits.push_back (best_bit);
      state = (state << 1) | best_bit;
    }

  /* remove termination */
  assert (decoded_bits.size() >= order);
  decoded_bits.resize (decoded_bits.size() - order);

  return decoded_bits;
}

int
main()
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

  vector<int> decoded_bits = conv_decode (coded_bits);
  printf ("output vector (k=%zd): ", decoded_bits.size());
  for (auto b : decoded_bits)
    printf ("%d", b);
  printf ("\n");
}
