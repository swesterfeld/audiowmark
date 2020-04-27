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

#include <string.h>
#include <stdio.h>

#include <array>
#include <regex>

#include "utils.hh"

using std::string;
using std::vector;
using std::regex;

class TSPacket
{
  std::array<unsigned char, 188> data;

public:
  bool
  read (FILE *file, Error& err)
  {
    err = Error::Code::NONE;

    size_t bytes_read = fread (data.data(), 1, data.size(), file);
    if (bytes_read == 0) /* probably eof */
      return false;

    if (bytes_read == data.size()) /* successful read */
      {
        if (data[0] == 'G')
          return true;

        err = Error ("bad packet sync while reading transport (.ts) packet");
        return false;
      }

    err = Error ("short read while reading transport stream (.ts) packet");
    return false;
  }
  Error
  write (FILE *file)
  {
    size_t bytes_written = fwrite (data.data(), 1, data.size(), file);
    if (bytes_written != data.size())
      return Error ("short write while writing transport stream (.ts) packet");

    return Error::Code::NONE;
  }
  enum Type { awmk_file, awmk_data };
  void
  clear (Type type)
  {
    std::fill (data.begin(), data.end(), 0);
    auto id = get_id_bytes (type);
    std::copy (id.begin(), id.end(), data.begin());
  }
  unsigned char&
  operator[] (size_t n)
  {
    return data[n];
  }
  bool
  has_id (Type type)
  {
    const auto idb = get_id_bytes (type);
    for (size_t i = 0; i < idb.size(); i++)
      {
        if (idb[i] != data[i])
          return false;
      }
    return true;
  }
  std::array<unsigned char, 12>
  get_id_bytes (Type type)
  {
    if (type == awmk_file)
     return { 'G', 0x1F, 0xFF, 0x10, 'A', 'W', 'M', 'K', 'f', 'i', 'l', 'e' };
    if (type == awmk_data)
      return { 'G', 0x1F, 0xFF, 0x10, 'A', 'W', 'M', 'K', 'd', 'a', 't', 'a' };
    return {0,};
  }
  constexpr size_t
  size()
  {
    return data.size(); // is constant
  }
};

Error
ts_append (const string& inname, const string& outname, const string& dataname)
{
  FILE *infile = fopen (inname.c_str(), "r");
  FILE *outfile = fopen (outname.c_str(), "w");

  while (!feof (infile))
    {
      TSPacket p;
      Error err;
      bool read_ok = p.read (infile, err);
      if (!read_ok)
        {
          if (err)
            return err;
        }
      else
        {
          err = p.write (outfile);
        }
    }
  vector<unsigned char> data;
  FILE *datafile = fopen (dataname.c_str(), "r");
  int c;
  while ((c = fgetc (datafile)) >= 0)
    data.push_back (c);

  char buf[1024];
  sprintf (buf, "%zd", data.size());
  string header = buf;
  header += ":" + dataname;
  header += '\0';
  for (size_t i = 0; i < header.size(); i++)
    data.insert (data.begin() + i, header[i]);

  TSPacket p_file;
  p_file.clear (TSPacket::awmk_file);
  size_t data_pos = 0;
  int pos = 12;
  while (data_pos < data.size())
    {
      p_file[pos++] = data[data_pos];
      if (pos == 188)
        {
          p_file.write (outfile);
          p_file.clear (TSPacket::awmk_data);
          pos = 12;
        }
      data_pos++;
    }

  if (pos != 12)
    {
      Error err = p_file.write (outfile);
    }

  return Error::Code::NONE;
}

class TSReader
{
public:
  struct Entry
  {
    string                filename;
    vector<unsigned char> data;
  };
private:
  struct Header
  {
    string filename;
    size_t data_size    = 0;
    size_t header_size  = 0;
  };
  vector<Entry> m_entries;
  bool parse_header (Header& header, const vector<unsigned char>& data);
public:
  Error load (const string& inname);
  const vector<Entry>& entries();
};

bool
TSReader::parse_header (Header& header, const vector<unsigned char>& data)
{
  for (size_t i = 0; i < data.size(); i++)
    {
      if (data[i] == 0) // header is terminated with one single 0 byte
        {
          string s = (const char *) (&data[0]);

          static const regex header_re ("([0-9]*):(.*)");
          std::smatch sm;
          if (regex_match (s, sm, header_re))
            {
              header.header_size = i + 1; // number of header bytes including 0 termination
              header.data_size = atoi (sm[1].str().c_str());
              header.filename = sm[2];
              return true;
            }
        }
    }
  return false;
}

Error
TSReader::load (const string& inname)
{
  FILE *infile = fopen (inname.c_str(), "r");

  vector<unsigned char> awmk_stream;
  Header header;
  bool header_valid = false;
  while (!feof (infile))
    {
      TSPacket p;
      Error err;
      bool read_ok = p.read (infile, err);
      if (!read_ok)
        {
          if (err)
            return err;
        }
      else
        {
          if (p.has_id (TSPacket::awmk_file) || p.has_id (TSPacket::awmk_data))
            {
              for (size_t i = 12; i < p.size(); i++)
                awmk_stream.push_back (p[i]);

              if (!header_valid)
                {
                  if (parse_header (header, awmk_stream))
                    header_valid = true;
                }
              // done? do we have enough bytes for the complete entry?
              if (header_valid && awmk_stream.size() >= (header.header_size + header.data_size))
                {
                  awmk_stream.erase (awmk_stream.begin(), awmk_stream.begin() + header.header_size);
                  awmk_stream.resize (header.data_size);

                  m_entries.push_back ({ header.filename, awmk_stream });

                  header_valid = false;
                  awmk_stream.clear();
                }
            }
        }
    }
  return Error::Code::NONE;
}

const vector<TSReader::Entry>&
TSReader::entries()
{
  return m_entries;
}

int
main (int argc, char **argv)
{
  if (argc == 5 && strcmp (argv[1], "append") == 0)
    {
      printf ("append: in=%s out=%s fn=%s\n", argv[2], argv[3], argv[4]);
      Error err = ts_append (argv[2], argv[3], argv[4]);
      if (err)
        {
          error ("ts_append: %s\n", err.message());
          return 1;
        }
    }
  else if (argc == 3 && strcmp (argv[1], "list") == 0)
    {
      TSReader reader;

      Error err = reader.load (argv[2]);
      for (auto entry : reader.entries())
        printf ("%s %zd\n", entry.filename.c_str(), entry.data.size());
    }
  else if (argc == 4 && strcmp (argv[1], "get") == 0)
    {
      TSReader reader;

      Error err = reader.load (argv[2]);
      for (auto entry : reader.entries())
        if (entry.filename == argv[3])
          fwrite (&entry.data[0], 1, entry.data.size(), stdout);
    }
  else
    {
      error ("testmpegts: error parsing command line arguments\n");
    }
}
