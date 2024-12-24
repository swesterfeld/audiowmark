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

#ifndef AUDIOWMARK_UTILS_HH
#define AUDIOWMARK_UTILS_HH

#include <vector>
#include <string>

#ifndef __STDC_FORMAT_MACROS
// some compilers/platforms (i.e. very old macOS) need this for macros like PRId64 (#61)
#define __STDC_FORMAT_MACROS
#endif
#include <cinttypes>

std::vector<int> bit_str_to_vec (const std::string& bits);
std::string      bit_vec_to_str (const std::vector<int>& bit_vec);

std::vector<unsigned char> hex_str_to_vec (const std::string& str);
std::string                vec_to_hex_str (const std::vector<unsigned char>& vec);

double get_time();

template<typename T>
inline const T&
bound (const T& min_value, const T& value, const T& max_value)
{
  return std::min (std::max (value, min_value), max_value);
}

// detect compiler
#if __clang__
  #define AUDIOWMARK_COMP_CLANG
  #define AUDIOWMARK_EXTRA_OPT
#elif __GNUC__ > 2
  #define AUDIOWMARK_COMP_GCC
  #define AUDIOWMARK_EXTRA_OPT __attribute__((optimize("-O3"))) /* enable auto vectorization, some functions benefit a lot from this */
#else
  #error "unsupported compiler"
#endif

#ifdef AUDIOWMARK_COMP_GCC
  #define AUDIOWMARK_PRINTF(format_idx, arg_idx)      __attribute__ ((__format__ (gnu_printf, format_idx, arg_idx)))
#else
  #define AUDIOWMARK_PRINTF(format_idx, arg_idx)      __attribute__ ((__format__ (__printf__, format_idx, arg_idx)))
#endif

/* bswap for g++ / clang++ - may need different implementation for other compilers */
static inline uint32_t
bswap32 (uint32_t i)
{
  return __builtin_bswap32 (i);
}

static inline uint64_t
bswap64 (uint64_t i)
{
  return __builtin_bswap64 (i);
}

void error (const char *format, ...) AUDIOWMARK_PRINTF (1, 2);
void warning (const char *format, ...) AUDIOWMARK_PRINTF (1, 2);
void info (const char *format, ...) AUDIOWMARK_PRINTF (1, 2);
void debug (const char *format, ...) AUDIOWMARK_PRINTF (1, 2);

enum class Log { ERROR = 3, WARNING = 2, INFO = 1, DEBUG = 0 };

void set_log_level (Log level);

std::string string_printf (const char *fmt, ...) AUDIOWMARK_PRINTF (1, 2);

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

class ScopedFile
{
  FILE *m_file;
public:
  ScopedFile (FILE *f) :
    m_file (f)
  {
  }
  ~ScopedFile()
  {
    if (m_file)
      fclose (m_file);
  }
};

#endif /* AUDIOWMARK_UTILS_HH */
