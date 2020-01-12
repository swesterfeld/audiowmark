#ifndef AUDIOWMARK_UTILS_HH
#define AUDIOWMARK_UTILS_HH

#include <vector>
#include <string>

std::vector<int> bit_str_to_vec (const std::string& bits);
std::string      bit_vec_to_str (const std::vector<int>& bit_vec);

std::vector<unsigned char> hex_str_to_vec (const std::string& str);
std::string                vec_to_hex_str (const std::vector<unsigned char>& vec);

template<typename T>
inline const T&
bound (const T& min_value, const T& value, const T& max_value)
{
  return std::min (std::max (value, min_value), max_value);
}

// detect compiler
#if __clang__
  #define AUDIOWMARK_COMP_CLANG
#elif __GNUC__ > 2
  #define AUDIOWMARK_COMP_GCC
#else
  #error "unsupported compiler"
#endif

#ifdef AUDIOWMARK_COMP_GCC
  #define AUDIOWMARK_PRINTF(format_idx, arg_idx)      __attribute__ ((__format__ (gnu_printf, format_idx, arg_idx)))
#else
  #define AUDIOWMARK_PRINTF(format_idx, arg_idx)      __attribute__ ((__format__ (__printf__, format_idx, arg_idx)))
#endif

void error (const char *format, ...) AUDIOWMARK_PRINTF (1, 2);
void warning (const char *format, ...) AUDIOWMARK_PRINTF (1, 2);
void info (const char *format, ...) AUDIOWMARK_PRINTF (1, 2);
void debug (const char *format, ...) AUDIOWMARK_PRINTF (1, 2);

enum class Log { ERROR = 3, WARNING = 2, INFO = 1, DEBUG = 0 };

void set_log_level (Log level);

class Error
{
public:
  enum class Code {
    NONE,
    STR
  };
  Error (Code code = Code::NONE) :
    m_code (code)
  {
    switch (code)
      {
        case Code::NONE:
          m_message = "OK";
          break;

        default:
          m_message = "Unknown error";
      }
  }
  explicit
  Error (const std::string& message) :
    m_code (Code::STR),
    m_message (message)
  {
  }
  Code
  code()
  {
    return m_code;
  }
  const char *
  message()
  {
    return m_message.c_str();
  }
  operator bool()
  {
    return m_code != Code::NONE;
  }
private:
  Code        m_code;
  std::string m_message;
};

#endif /* AUDIOWMARK_UTILS_HH */