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

#include "utils.hh"
#include "convcode.hh"

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

constexpr  auto         ab_generators = std::array<unsigned,12>
  {
    066561, 075211, 071545, 054435, 063635, 052475,
    063543, 075307, 052547, 045627, 067657, 051757
  };

constexpr  unsigned int ab_rate       = ab_generators.size();
constexpr  unsigned int order         = 15;

/*
constexpr  unsigned int order       = 9;
constexpr  auto         generators  = std::array<unsigned,3> { 0557, 0663, 0711 };

constexpr  unsigned int order       = 14;
constexpr  auto         generators  = std::array<unsigned,3> { 021645, 035661, 037133 };

constexpr  unsigned int order       = 18;
constexpr  auto         generators  = std::array<unsigned,3> { 0552137, 0614671, 0772233 };
*/

constexpr  unsigned int state_count = (1 << order);
constexpr  unsigned int state_mask  = (1 << order) - 1;

size_t
conv_code_size (ConvBlockType block_type, size_t msg_size)
{
  switch (block_type)
    {
      case ConvBlockType::a:
      case ConvBlockType::b:  return (msg_size + order) * ab_rate / 2;
      case ConvBlockType::ab: return (msg_size + order) * ab_rate;
      default:                assert (false);
    }
}

vector<unsigned>
get_block_type_generators (ConvBlockType block_type)
{
  vector<unsigned> generators;

  if (block_type == ConvBlockType::a)
    {
      for (unsigned int i = 0; i < ab_rate / 2; i++)
        generators.push_back (ab_generators[i * 2]);
    }
  else if (block_type == ConvBlockType::b)
    {
      for (unsigned int i = 0; i < ab_rate / 2; i++)
        generators.push_back (ab_generators[i * 2 + 1]);
    }
  else
    {
      assert (block_type == ConvBlockType::ab);
      generators.assign (ab_generators.begin(), ab_generators.end());
    }
  return generators;
}

vector<int>
conv_encode (ConvBlockType block_type, const vector<int>& in_bits)
{
  auto generators = get_block_type_generators (block_type);

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

/* decode using viterbi algorithm */
vector<int>
conv_decode_soft (ConvBlockType block_type, const vector<float>& coded_bits, float *error_out)
{
  auto generators = get_block_type_generators (block_type);
  unsigned int rate = generators.size();
  vector<int> decoded_bits;

  assert (coded_bits.size() % rate == 0);

  struct StateEntry
  {
    int   last_state;
    float delta;
    int   bit;
  };
  vector<vector<StateEntry>> error_count;
  for (size_t i = 0; i < coded_bits.size() + rate; i += rate) /* 1 extra element */
    error_count.emplace_back (state_count, StateEntry {0, -1, 0});

  error_count[0][0].delta = 0; /* start state */

  /* precompute state -> output bits table */
  vector<float> state2bits;
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
                  int   new_state = ((state << 1) | bit) & state_mask;

                  float delta = old_table[state].delta;
                  int   sbit_pos = new_state * rate;

                  for (size_t p = 0; p < rate; p++)
                    {
                      const float cbit = coded_bits[i + p];
                      const float sbit = state2bits[sbit_pos + p];

                      /* decoding error weight for this bit; if input is only 0.0 and 1.0, this is the hamming distance */
                      delta += (cbit - sbit) * (cbit - sbit);
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

  unsigned int state = 0;
  if (error_out)
    *error_out = error_count.back()[state].delta / coded_bits.size();
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

vector<int>
conv_decode_hard (ConvBlockType block_type, const vector<int>& coded_bits)
{
  /* for the final application, we always want soft decoding, so we don't
   * special case hard decoding here, so this will be a little slower than
   * necessary
   */
  vector<float> soft_bits;
  for (auto b : coded_bits)
    soft_bits.push_back (b ? 1.0f : 0.0f);

  return conv_decode_soft (block_type, soft_bits);
}

void
conv_print_table (ConvBlockType block_type)
{
  vector<int> bits (100);
  bits[0] = 1;

  vector<int> out_bits = conv_encode (block_type, bits);

  auto generators = get_block_type_generators (block_type);
  unsigned int rate = generators.size();

  for (unsigned int r = 0; r < rate; r++)
    {
      for (unsigned int i = 0; i < order; i++)
        printf ("%s%d", i == 0 ? "" : " ", out_bits[i * rate + r]);
      printf ("\n");
    }
}
