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

#ifndef AUDIOWMARK_SHORT_CODE_HH
#define AUDIOWMARK_SHORT_CODE_HH

#include <vector>
#include <string>

#include "convcode.hh"

size_t           code_size (ConvBlockType block_type, size_t msg_size);
std::vector<int> code_encode (ConvBlockType block_type, const std::vector<int>& in_bits);
std::vector<int> code_decode_soft (ConvBlockType block_type, const std::vector<float>& coded_bits, float *error_out = nullptr);

size_t           short_code_size (ConvBlockType block_type, size_t msg_size);
std::vector<int> short_encode (ConvBlockType block_type, const std::vector<int>& in_bits);
std::vector<int> short_decode_soft (ConvBlockType block_type, const std::vector<float>& coded_bits, float *error_out = nullptr);

std::vector<int> short_encode_blk (const std::vector<int>& in_bits);
std::vector<int> short_decode_blk (const std::vector<int>& coded_bits);
size_t           short_code_init (size_t k);

#endif /* AUDIOWMARK_SHORT_CODE_HH */
