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

  virtual void to_raw   (const std::vector<float>& samples, std::vector<unsigned char>& bytes) = 0;
  virtual void from_raw (const std::vector<unsigned char>& bytes, std::vector<float>& samples) = 0;
};

#endif /* AUDIOWMARK_RAW_CONVERTER_HH */
