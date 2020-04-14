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
