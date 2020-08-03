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

#ifndef AUDIOWMARK_LIMITER_HH
#define AUDIOWMARK_LIMITER_HH

#include <vector>
#include <sys/types.h>

class Limiter
{
  float ceiling           = 1;
  float block_max_last    = 0;
  float block_max_current = 0;
  float block_max_next    = 0;
  uint  block_size        = 0;
  uint  n_channels        = 0;
  uint  sample_rate       = 0;

  std::vector<float> buffer;
  void process_block (const float *in, float *out);
  float block_max (const float *in);
  void debug_scale (float scale);
public:
  Limiter (int n_channels, int sample_rate);

  void set_block_size_ms (int value_ms);
  void set_ceiling (float ceiling);

  std::vector<float> process (const std::vector<float>& samples);
  size_t             skip (size_t zeros);
  std::vector<float> flush();
};

#endif /* AUDIOWMARK_LIMITER_HH */
