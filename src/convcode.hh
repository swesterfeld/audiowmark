#ifndef AUDIOWMARK_CONV_CODE_HH
#define AUDIOWMARK_CONV_CODE_HH

#include <vector>
#include <string>

enum class ConvBlockType { a, b, ab };

size_t           conv_code_size (ConvBlockType block_type, size_t msg_size);
std::vector<int> conv_encode (ConvBlockType block_type, const std::vector<int>& in_bits);
std::vector<int> conv_decode_hard (ConvBlockType block_type, const std::vector<int>& coded_bits);
std::vector<int> conv_decode_soft (ConvBlockType block_type, const std::vector<float>& coded_bits, float *error_out = nullptr);

void             conv_print_table (ConvBlockType block_type);

#endif /* AUDIOWMARK_CONV_CODE_HH */
