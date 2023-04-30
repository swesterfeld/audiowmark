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

#include "utils.hh"
#include "stdarg.h"

#include <sys/time.h>

using std::vector;
using std::string;

double
get_time()
{
  /* return timestamp in seconds as double */
  timeval tv;
  gettimeofday (&tv, 0);

  return tv.tv_sec + tv.tv_usec / 1000000.0;
}

static unsigned char
from_hex_nibble (char c)
{
  int uc = (unsigned char)c;

  if (uc >= '0' && uc <= '9') return uc - (unsigned char)'0';
  if (uc >= 'a' && uc <= 'f') return uc + 10 - (unsigned char)'a';
  if (uc >= 'A' && uc <= 'F') return uc + 10 - (unsigned char)'A';

  return 16;	// error
}

vector<int>
bit_str_to_vec (const string& bits)
{
  vector<int> bitvec;
  for (auto nibble : bits)
    {
      unsigned char c = from_hex_nibble (nibble);
      if (c >= 16)
        return vector<int>(); // error

      bitvec.push_back ((c & 8) > 0);
      bitvec.push_back ((c & 4) > 0);
      bitvec.push_back ((c & 2) > 0);
      bitvec.push_back ((c & 1) > 0);
    }
  return bitvec;
}

string
bit_vec_to_str (const vector<int>& bit_vec)
{
  string bit_str;

  for (size_t pos = 0; pos + 3 < bit_vec.size(); pos += 4) // convert only groups of 4 bits
    {
      int nibble = 0;
      for (int j = 0; j < 4; j++)
        {
          if (bit_vec[pos + j])
            {
              // j == 0 has the highest value, then 1, 2, 3 (lowest)
              nibble |= 1 << (3 - j);
            }
        }
      const char *to_hex = "0123456789abcdef";
      bit_str += to_hex[nibble];
    }
  return bit_str;
}

vector<unsigned char>
hex_str_to_vec (const string& str)
{
  vector<unsigned char> result;

  if ((str.size() % 2) != 0) // even length
    return vector<unsigned char>();

  for (size_t i = 0; i < str.size() / 2; i++)
    {
      unsigned char h = from_hex_nibble (str[i * 2]);
      unsigned char l = from_hex_nibble (str[i * 2 + 1]);
      if (h >= 16 || l >= 16)
        return vector<unsigned char>();

      result.push_back ((h << 4) + l);
    }

  return result;
}

string
vec_to_hex_str (const vector<unsigned char>& vec)
{
  string s;
  for (auto byte : vec)
    s += string_printf ("%02x", byte);

  return s;
}

static string
string_vprintf (const char *format, va_list vargs)
{
  string s;

  char *str = NULL;
  if (vasprintf (&str, format, vargs) >= 0 && str)
    {
      s = str;
      free (str);
    }
  else
    s = format;

  return s;
}

string
string_printf (const char *format, ...)
{
  va_list ap;

  va_start (ap, format);
  string s = string_vprintf (format, ap);
  va_end (ap);

  return s;
}

static Log log_level = Log::INFO;

void
set_log_level (Log level)
{
  log_level = level;
}

static void
logv (Log log, const char *format, va_list vargs)
{
  if (log >= log_level)
    {
      string s = string_vprintf (format, vargs);

      /* could support custom log function here */
      fprintf (stderr, "%s", s.c_str());
      fflush (stderr);
    }
}

void
error (const char *format, ...)
{
  va_list ap;

  va_start (ap, format);
  logv (Log::ERROR, format, ap);
  va_end (ap);
}

void
warning (const char *format, ...)
{
  va_list ap;

  va_start (ap, format);
  logv (Log::WARNING, format, ap);
  va_end (ap);
}

void
info (const char *format, ...)
{
  va_list ap;

  va_start (ap, format);
  logv (Log::INFO, format, ap);
  va_end (ap);
}

void
debug (const char *format, ...)
{
  va_list ap;

  va_start (ap, format);
  logv (Log::DEBUG, format, ap);
  va_end (ap);
}
