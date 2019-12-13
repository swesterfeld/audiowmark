#ifndef AUDIOWMARK_RAW_CONVERTER_HH
#define AUDIOWMARK_RAW_CONVERTER_HH

#include "rawinputstream.hh"

class RawConverter
{
public:
  static RawConverter *create (const RawFormat& raw_format, Error& error);

  virtual void to_raw   (const std::vector<float>& samples, std::vector<unsigned char>& bytes) = 0;
  virtual void from_raw (const std::vector<unsigned char>& bytes, std::vector<float>& samples) = 0;
};

#endif /* AUDIOWMARK_RAW_CONVERTER_HH */
