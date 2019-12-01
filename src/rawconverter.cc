#include "rawconverter.hh"

#include <math.h>

using std::vector;

template<int BIT_DEPTH>
class RawConverterImpl : public RawConverter
{
public:
  void to_raw (const std::vector<float>& samples, std::vector<unsigned char>& bytes);
  void from_raw (const std::vector<unsigned char>& bytes, std::vector<float>& samples);
};

RawConverter *
RawConverter::create (const RawFormat& raw_format, Error& error)
{
  error = Error::Code::NONE;
  switch (raw_format.bit_depth())
    {
      case 16: return new RawConverterImpl<16>();
      case 24: return new RawConverterImpl<24>();
      default: error = Error ("unsupported bit depth");
               return nullptr;
    }
}

template<int BIT_DEPTH>
void
RawConverterImpl<BIT_DEPTH>::to_raw (const vector<float>& samples, vector<unsigned char>& output_bytes)
{
  output_bytes.resize (BIT_DEPTH / 8 * samples.size());

  for (size_t i = 0; i < samples.size(); i++)
    {
      const double norm      =  0x80000000LL;
      const double min_value = -0x80000000LL;
      const double max_value =  0x7FFFFFFF;

      const int    sample = lrint (bound<double> (min_value, samples[i] * norm, max_value));

      if (BIT_DEPTH == 16)
        {
          // write 16-bit little endian value
          output_bytes[i * 2]     = sample >> 16;
          output_bytes[i * 2 + 1] = sample >> 24;
        }
      else if (BIT_DEPTH == 24)
        {
          // write 24-bit little endian value
          output_bytes[i * 3]     = sample >> 8;
          output_bytes[i * 3 + 1] = sample >> 16;
          output_bytes[i * 3 + 2] = sample >> 24;
        }
    }
}

template<int BIT_DEPTH>
void
RawConverterImpl<BIT_DEPTH>::from_raw (const vector<unsigned char>& input_bytes, vector<float>& samples)
{
  const unsigned char *ptr = input_bytes.data();

  samples.resize (input_bytes.size() / (BIT_DEPTH / 8));
  const double norm = 1.0 / 0x80000000LL;
  for (size_t i = 0; i < samples.size(); i++)
    {
      int s32 = (ptr[1] << 24) + (ptr[0] << 16);
      samples[i] = s32 * norm;
      ptr += 2;
    }
}
