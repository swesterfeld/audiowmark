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
#include "mpegts.hh"

using std::string;
using std::vector;
using std::map;
using std::regex;

int
main (int argc, char **argv)
{
  if (argc == 5 && strcmp (argv[1], "append") == 0)
    {
      printf ("append: in=%s out=%s fn=%s\n", argv[2], argv[3], argv[4]);
      TSWriter writer;

      writer.append_file (argv[4], argv[4]);
      Error err = writer.process (argv[2], argv[3]);
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
  else if (argc == 3 && strcmp (argv[1], "vars") == 0)
    {
      TSReader reader;

      Error err = reader.load (argv[2]);
      map<string, string> vars = reader.parse_vars ("vars");
      for (auto v : vars)
        printf ("%s=%s\n", v.first.c_str(), v.second.c_str());
    }
  else if (argc == 3 && strcmp (argv[1], "perf") == 0)
    {
      for (int i = 0; i < 1000; i++)
        {
          TSReader reader;

          Error err = reader.load (argv[2]);
          if (i == 42)
            for (auto entry : reader.entries())
              printf ("%s %zd\n", entry.filename.c_str(), entry.data.size());
        }
    }
  else
    {
      error ("testmpegts: error parsing command line arguments\n");
    }
}
