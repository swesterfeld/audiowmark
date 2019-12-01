#include "rawconverter.hh"

#include <math.h>

using std::vector;

class RawConverterImpl : public RawConverter
{
public:
  void to_raw (const std::vector<float>& samples, std::vector<unsigned char>& bytes);
  void
  from_raw (const std::vector<unsigned char>& bytes, std::vector<float>& samples)
  {
  }
};

RawConverter *
RawConverter::create (const RawFormat& raw_format, Error& error)
{
  error = Error::Code::NONE;
  return new RawConverterImpl;
}

void
RawConverterImpl::to_raw (const vector<float>& samples, vector<unsigned char>& output_bytes)
{
  output_bytes.resize (2 * samples.size());

  for (size_t i = 0; i < samples.size(); i++)
    {
      const double norm      =  0x80000000LL;
      const double min_value = -0x80000000LL;
      const double max_value =  0x7FFFFFFF;

      const int    sample = lrint (bound<double> (min_value, samples[i] * norm, max_value));

      // write 16-bit little endian value
      output_bytes[i * 2]     = sample >> 16;
      output_bytes[i * 2 + 1] = sample >> 24;
    }
}
