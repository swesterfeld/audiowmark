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

#ifndef AUDIOWMARK_RAW_CONVERTER_HH
#define AUDIOWMARK_RAW_CONVERTER_HH

#include "rawinputstream.hh"

class RawConverter
{
public:
  static RawConverter *create (const RawFormat& raw_format, Error& error);

  virtual ~RawConverter() = 0;

  virtual void to_raw   (const float *samples, unsigned char *bytes, size_t n_samples) = 0;
  virtual void from_raw (const unsigned char *bytes, float *samples, size_t n_samples) = 0;
};

template<int BITS>
static inline int
float_to_int_clip (float f)
{
  const int64_t inorm = (1LL << (BITS - 1));
  const float min_value = -inorm;
  const float max_value =  inorm - 1;
  const float norm      =  inorm;
  const float snorm     = f * norm;

  if (snorm >= max_value)
    return inorm - 1;
  else if (snorm <= min_value)
    return -inorm;
  else
    return snorm;
}

static inline float
float_clip (float f)
{
  const float min_value = -1;
  const float max_value =  1;

  if (f >= max_value)
    return max_value;
  else if (f <= min_value)
    return min_value;
  else
    return f;
}

#endif /* AUDIOWMARK_RAW_CONVERTER_HH */
