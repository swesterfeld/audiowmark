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

#ifndef AUDIOWMARK_MPEGTS_HH
#define AUDIOWMARK_MPEGTS_HH

class TSReader
{
public:
  struct Entry
  {
    std::string                filename;
    std::vector<unsigned char> data;
  };
private:
  struct Header
  {
    std::string filename;
    size_t      data_size = 0;
  };
  std::vector<Entry> m_entries;
  bool parse_header (Header& header, std::vector<unsigned char>& data);
public:
  Error load (const std::string& inname);
  const std::vector<Entry>& entries();
};

Error ts_append (const std::string& inname, const std::string& outname, const std::string& dataname);

#endif /* AUDIOWMARK_MPEGTS_HH */
