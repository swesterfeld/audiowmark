#ifndef AUDIOWMARK_UTILS_HH
#define AUDIOWMARK_UTILS_HH

#include <vector>
#include <string>

std::vector<int> bit_str_to_vec (const std::string& bits);
std::string      bit_vec_to_str (const std::vector<int>& bit_vec);

std::vector<unsigned char> hex_str_to_vec (const std::string& str);
std::string                vec_to_hex_str (const std::vector<unsigned char>& vec);

#endif /* AUDIOWMARK_UTILS_HH */
