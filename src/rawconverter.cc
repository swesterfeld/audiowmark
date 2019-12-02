#include "rawconverter.hh"

#include <math.h>

using std::vector;

template<int BIT_DEPTH, RawFormat::Endian ENDIAN>
class RawConverterImpl : public RawConverter
{
public:
  void to_raw (const std::vector<float>& samples, std::vector<unsigned char>& bytes);
  void from_raw (const std::vector<unsigned char>& bytes, std::vector<float>& samples);
};

template<int BIT_DEPTH>
static RawConverter *
create_with_bits (const RawFormat& raw_format, Error& error)
{
  switch (raw_format.endian())
    {
      case RawFormat::LITTLE: return new RawConverterImpl<BIT_DEPTH, RawFormat::LITTLE>();
      case RawFormat::BIG:    return new RawConverterImpl<BIT_DEPTH, RawFormat::BIG>();
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
      default: error = Error ("unsupported bit depth");
               return nullptr;
    }
}

template<int BIT_DEPTH, RawFormat::Endian ENDIAN>
constexpr std::array<int, 3>
make_endian_shift ()
{
  if (BIT_DEPTH == 16)
    {
      if (ENDIAN == RawFormat::Endian::LITTLE)
        return { 16, 24, -1 };
      else
        return { 24, 16, -1 };
    }
  if (BIT_DEPTH == 24)
    {
      if (ENDIAN == RawFormat::Endian::LITTLE)
        return {  8, 16, 24 };
      else
        return { 24, 16, 8 };
    }
}

template<int BIT_DEPTH, RawFormat::Endian ENDIAN>
void
RawConverterImpl<BIT_DEPTH, ENDIAN>::to_raw (const vector<float>& samples, vector<unsigned char>& output_bytes)
{
  constexpr int  sample_width = BIT_DEPTH / 8;
  constexpr auto eshift = make_endian_shift<BIT_DEPTH, ENDIAN>();

  output_bytes.resize (sample_width * samples.size());

  unsigned char *ptr = output_bytes.data();

  for (size_t i = 0; i < samples.size(); i++)
    {
      const double norm      =  0x80000000LL;
      const double min_value = -0x80000000LL;
      const double max_value =  0x7FFFFFFF;

      const int    sample = lrint (bound<double> (min_value, samples[i] * norm, max_value));

      if (eshift[0] >= 0)
        ptr[0] = sample >> eshift[0];
      if (eshift[1] >= 0)
        ptr[1] = sample >> eshift[1];
      if (eshift[2] >= 0)
        ptr[2] = sample >> eshift[2];

      ptr += sample_width;
    }
}

template<int BIT_DEPTH, RawFormat::Endian ENDIAN>
void
RawConverterImpl<BIT_DEPTH, ENDIAN>::from_raw (const vector<unsigned char>& input_bytes, vector<float>& samples)
{
  const unsigned char *ptr = input_bytes.data();
  constexpr int sample_width = BIT_DEPTH / 8;
  constexpr auto eshift = make_endian_shift<BIT_DEPTH, ENDIAN>();

  samples.resize (input_bytes.size() / sample_width);
  const double norm = 1.0 / 0x80000000LL;
  for (size_t i = 0; i < samples.size(); i++)
    {
      int s32 = 0;

      if (eshift[0] >= 0)
        s32 += ptr[0] << eshift[0];
      if (eshift[1] >= 0)
        s32 += ptr[1] << eshift[1];
      if (eshift[2] >= 0)
        s32 += ptr[2] << eshift[2];

      samples[i] = s32 * norm;
      ptr += sample_width;
    }
}
