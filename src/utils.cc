#include "utils.hh"

using std::vector;
using std::string;

static unsigned char
from_hex_nibble (char c)
{
  int uc = (unsigned char)c;

  if (uc >= '0' && uc <= '9') return uc - (unsigned char)'0';
  if (uc >= 'a' && uc <= 'f') return uc + 10 - (unsigned char)'a';
  if (uc >= 'A' && uc <= 'F') return uc + 10 - (unsigned char)'A';

  return 16;	// error
}

vector<int>
bit_str_to_vec (const string& bits)
{
  vector<int> bitvec;
  for (auto nibble : bits)
    {
      unsigned char c = from_hex_nibble (nibble);
      if (c >= 16)
        return vector<int>(); // error

      bitvec.push_back ((c & 8) > 0);
      bitvec.push_back ((c & 4) > 0);
      bitvec.push_back ((c & 2) > 0);
      bitvec.push_back ((c & 1) > 0);
    }
  return bitvec;
}

string
bit_vec_to_str (const vector<int>& bit_vec)
{
  string bit_str;

  for (size_t pos = 0; pos + 3 < bit_vec.size(); pos += 4) // convert only groups of 4 bits
    {
      int nibble = 0;
      for (int j = 0; j < 4; j++)
        {
          if (bit_vec[pos + j])
            {
              // j == 0 has the highest value, then 1, 2, 3 (lowest)
              nibble |= 1 << (3 - j);
            }
        }
      const char *to_hex = "0123456789abcdef";
      bit_str += to_hex[nibble];
    }
  return bit_str;
}

