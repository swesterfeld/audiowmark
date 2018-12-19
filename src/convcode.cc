#include "utils.hh"

#include <array>
#include <algorithm>
#include <assert.h>

using std::vector;
using std::string;
using std::min;

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
/*
constexpr  unsigned int order       = 14;
constexpr  auto         generators  = std::array<unsigned,3> { 021645, 035661, 037133 };
*/
/*
constexpr  unsigned int order       = 18;
constexpr  auto         generators  = std::array<unsigned,3> { 0552137, 0614671, 0772233 };
*/

constexpr  unsigned int state_count = (1 << order);
constexpr  unsigned int state_mask  = (1 << order) - 1;

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

  struct StateEntry
  {
    int last_state;
    int delta;
    int bit;
  };
  vector<vector<StateEntry>> error_count;
  for (size_t i = 0; i < coded_bits.size() + rate; i += rate) /* 1 extra element */
    {
      vector<StateEntry> state_table;
      for (unsigned s = 0; s < state_count; s++)
        {
          if (s == 0 && i == 0)
            state_table.push_back ({0,0});
          else
            state_table.push_back ({-1,-1});
        }
      error_count.push_back (state_table);
    }

  for (size_t i = 0; i < coded_bits.size(); i += rate)
    {
      vector<StateEntry>& old_table = error_count[i / rate];
      vector<StateEntry>& new_table = error_count[i / rate + 1];

      for (unsigned int state = 0; state < state_count; state++)
        {
          /* this check enforces that we only consider states reachable from state=0 at time=0*/
          if (old_table[state].delta >= 0)
            {
              for (int bit = 0; bit < 2; bit++)
                {
                  int new_state = ((state << 1) | bit) & state_mask;

                  int delta = old_table[state].delta;
                  for (size_t p = 0; p < generators.size(); p++)
                    {
                      int out_bit = parity (new_state & generators[p]);
                      if (out_bit != coded_bits[i + p])
                        delta += 1;
                    }
                  if (delta < new_table[new_state].delta || new_table[new_state].delta < 0) /* better match with this link? replace entry */
                    {
                      new_table[new_state].delta      = delta;
                      new_table[new_state].last_state = state;
                      new_table[new_state].bit        = bit;
                    }
                }
            }
        }
    }

  size_t idx = error_count.size() - 1;
  unsigned int state = 0;
  do
    {
      decoded_bits.push_back (error_count[idx][state].bit);

      state = error_count[idx][state].last_state;
      idx--;
    }
  while (idx > 0);
  std::reverse (decoded_bits.begin(), decoded_bits.end());

  /* remove termination */
  assert (decoded_bits.size() >= order);
  decoded_bits.resize (decoded_bits.size() - order);

  return decoded_bits;
}

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

      vector<int> decoded_bits = conv_decode (coded_bits);
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
          constexpr int test_size = 1000;

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

              vector<int> decoded_bits = conv_decode (coded_bits);

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
}
