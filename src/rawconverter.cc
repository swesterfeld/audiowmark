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

#include "rawconverter.hh"
#include "config.h"

#include <array>

#include <math.h>
#include <assert.h>

using std::vector;

RawConverter::~RawConverter()
{
}

template<int BIT_DEPTH, RawFormat::Endian ENDIAN, RawFormat::Encoding ENCODING>
class RawConverterImpl final : public RawConverter
{
public:
  void to_raw (const float *samples, unsigned char *bytes, size_t n_samples) override AUDIOWMARK_EXTRA_OPT;
  void from_raw (const unsigned char *bytes, float *samples, size_t n_samples) override AUDIOWMARK_EXTRA_OPT;
};

template<int BIT_DEPTH, RawFormat::Endian ENDIAN>
static RawConverter *
create_with_bits_endian (const RawFormat& raw_format, Error& error)
{
  switch (raw_format.encoding())
    {
      case RawFormat::SIGNED:   return new RawConverterImpl<BIT_DEPTH, ENDIAN, RawFormat::SIGNED>();
      case RawFormat::UNSIGNED: return new RawConverterImpl<BIT_DEPTH, ENDIAN, RawFormat::UNSIGNED>();
    }
  error = Error ("unsupported encoding");
  return nullptr;
}

template<int BIT_DEPTH>
static RawConverter *
create_with_bits (const RawFormat& raw_format, Error& error)
{
  switch (raw_format.endian())
    {
      case RawFormat::LITTLE: return create_with_bits_endian <BIT_DEPTH, RawFormat::LITTLE> (raw_format, error);
      case RawFormat::BIG:    return create_with_bits_endian <BIT_DEPTH, RawFormat::BIG> (raw_format, error);
    }
  error = Error ("unsupported endianness");
  return nullptr;
}

RawConverter *
RawConverter::create (const RawFormat& raw_format, Error& error)
{
  error = Error::Code::NONE;
  switch (raw_format.bit_depth())
    {
      case 16: return create_with_bits<16> (raw_format, error);
      case 24: return create_with_bits<24> (raw_format, error);
      case 32: return create_with_bits<32> (raw_format, error);
      default: error = Error ("unsupported bit depth");
               return nullptr;
    }
}

template<int BIT_DEPTH, RawFormat::Endian ENDIAN>
constexpr std::array<int, 4>
make_endian_shift ()
{
  if (BIT_DEPTH == 16)
    {
      if (ENDIAN == RawFormat::Endian::LITTLE)
        return { 16, 24, -1, -1 };
      else
        return { 24, 16, -1, -1 };
    }
  if (BIT_DEPTH == 24)
    {
      if (ENDIAN == RawFormat::Endian::LITTLE)
        return {  8, 16, 24, -1 };
      else
        return { 24, 16, 8, -1 };
    }
  if (BIT_DEPTH == 32)
    {
      if (ENDIAN == RawFormat::Endian::LITTLE)
        return {  0, 8, 16, 24 };
      else
        return { 24, 16, 8, 0 };
    }
}

template<int BIT_DEPTH, RawFormat::Endian ENDIAN, RawFormat::Encoding ENCODING>
void
RawConverterImpl<BIT_DEPTH, ENDIAN, ENCODING>::to_raw (const float *samples, unsigned char *output_bytes, size_t n_samples)
{
  constexpr int  sample_width = BIT_DEPTH / 8;
  constexpr auto eshift = make_endian_shift<BIT_DEPTH, ENDIAN>();
  constexpr unsigned char sign_flip = ENCODING == RawFormat::SIGNED ? 0x00 : 0x80;
#ifdef WORDS_BIGENDDIAN
  constexpr bool native_endian = ENDIAN == RawFormat::BIG;
#else
  constexpr bool native_endian = ENDIAN == RawFormat::LITTLE;
#endif
  assert ((uintptr_t (output_bytes) & 3) == 0); // ensure alignment for optimized 32 bit native endian version

  unsigned char *ptr = output_bytes;

  for (size_t i = 0; i < n_samples; i++)
    {
      if (native_endian && ENCODING == RawFormat::SIGNED && BIT_DEPTH == 32)
        ((int32_t *)ptr)[i] = float_to_int_clip<32> (samples[i]);
      else if (native_endian && ENCODING == RawFormat::SIGNED && BIT_DEPTH == 16)
        ((int16_t *)ptr)[i] = float_to_int_clip<16> (samples[i]);
      else
        {
          int sample = float_to_int_clip<32> (samples[i]);

          if (eshift[0] >= 0)
            ptr[0] = (sample >> eshift[0]) ^ sign_flip;
          if (eshift[1] >= 0)
            ptr[1] = sample >> eshift[1];
          if (eshift[2] >= 0)
            ptr[2] = sample >> eshift[2];
          if (eshift[3] >= 0)
            ptr[3] = sample >> eshift[3];

          ptr += sample_width;
        }
    }
}

template<int BIT_DEPTH, RawFormat::Endian ENDIAN, RawFormat::Encoding ENCODING>
void
RawConverterImpl<BIT_DEPTH, ENDIAN, ENCODING>::from_raw (const unsigned char *input_bytes, float *samples, size_t n_samples)
{
  const unsigned char *ptr = input_bytes;
  constexpr int sample_width = BIT_DEPTH / 8;
  constexpr auto eshift = make_endian_shift<BIT_DEPTH, ENDIAN>();
  constexpr unsigned char sign_flip = ENCODING == RawFormat::SIGNED ? 0x00 : 0x80;
#ifdef WORDS_BIGENDDIAN
  constexpr bool native_endian = ENDIAN == RawFormat::BIG;
#else
  constexpr bool native_endian = ENDIAN == RawFormat::LITTLE;
#endif
  assert ((uintptr_t (input_bytes) & 3) == 0); // ensure alignment for optimized 32 bit native endian version

  const float norm = 1.0 / 0x80000000LL;
  for (size_t i = 0; i < n_samples; i++)
    {
      if (native_endian && ENCODING == RawFormat::SIGNED && BIT_DEPTH == 32)
        samples[i] = ((int32_t *)ptr)[i] * norm;
      else if (native_endian && ENCODING == RawFormat::SIGNED && BIT_DEPTH == 16)
        samples[i] = ((int16_t *)ptr)[i] * float (1.0 / 0x8000);
      else
        {
          int s32 = 0;

          if (eshift[0] >= 0)
            s32 += (ptr[0] ^ sign_flip) << eshift[0];
          if (eshift[1] >= 0)
            s32 += ptr[1] << eshift[1];
          if (eshift[2] >= 0)
            s32 += ptr[2] << eshift[2];
          if (eshift[3] >= 0)
            s32 += ptr[3] << eshift[3];
          samples[i] = s32 * norm;

          ptr += sample_width;
        }

    }
}
