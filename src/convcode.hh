#ifndef AUDIOWMARK_CONV_CODE_HH
#define AUDIOWMARK_CONV_CODE_HH

#include <vector>
#include <string>

std::vector<int> conv_encode (const std::vector<int>& in_bits);
std::vector<int> conv_decode (const std::vector<int>& coded_bits);

#endif /* AUDIOWMARK_CONV_CODE_HH */
