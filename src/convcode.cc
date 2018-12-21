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
/*
constexpr  unsigned int order       = 9;
constexpr  auto         generators  = std::array<unsigned,3> { 0557, 0663, 0711 };
*/
constexpr  unsigned int order       = 14;
constexpr  auto         generators  = std::array<unsigned,3> { 021645, 035661, 037133 };
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
    error_count.emplace_back (state_count, StateEntry {0, -1, 0});

  error_count[0][0].delta = 0; /* start state */

  /* precompute state -> output bits table */
  vector<uint8_t> state2bits;
  for (unsigned int state = 0; state < state_count; state++)
    {
      for (size_t p = 0; p < generators.size(); p++)
        {
          int out_bit = parity (state & generators[p]);
          state2bits.push_back (out_bit);
        }
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
                  int sbit_pos = new_state * rate;

                  /* hamming distance between produced bits and coded bits */
                  for (size_t p = 0; p < generators.size(); p++)
                    delta += state2bits[sbit_pos + p] ^ coded_bits[i + p];

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

  unsigned int state = 0;
  for (size_t idx = error_count.size() - 1; idx > 0; idx--)
    {
      decoded_bits.push_back (error_count[idx][state].bit);

      state = error_count[idx][state].last_state;
    }
  std::reverse (decoded_bits.begin(), decoded_bits.end());

  /* remove termination */
  assert (decoded_bits.size() >= order);
  decoded_bits.resize (decoded_bits.size() - order);

  return decoded_bits;
}
