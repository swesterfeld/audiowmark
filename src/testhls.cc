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

#include <regex>

#include "utils.hh"
#include "mpegts.hh"
#include "wavdata.hh"

using std::string;
using std::regex;
using std::vector;
using std::map;

Error
ff_decode (const string& filename, WavData& out_wav_data)
{
  FILE *tmp_file = tmpfile();
  ScopedFile tmp_file_s (tmp_file);
  string tmp_file_name = string_printf ("/dev/fd/%d", fileno (tmp_file));

  string cmd = string_printf ("ffmpeg -v error -y -i '%s' -f wav %s", filename.c_str(), tmp_file_name.c_str());
  system (cmd.c_str());

  Error err = out_wav_data.load (tmp_file_name);
  return err;
}

Error
ff_encode (const WavData& wav_data, const string& filename, size_t start_pos)
{
  FILE *tmp_file = tmpfile();
  ScopedFile tmp_file_s (tmp_file);
  string tmp_file_name = string_printf ("/dev/fd/%d", fileno (tmp_file));

  Error err = wav_data.save (tmp_file_name);

  string cmd = string_printf ("ffmpeg -v error -y -i %s -f mpegts -af asetpts='(%zd+N)/SR/TB' -c:a aac '%s'",
                              tmp_file_name.c_str(), start_pos, filename.c_str());
  system (cmd.c_str());

  return Error::Code::NONE;
}

int
hls_embed_context (const string& in_dir, const string& out_dir, const string& filename)
{
  string in_name = in_dir + "/" + filename;
  FILE *in_file = fopen (in_name.c_str(), "r");
  ScopedFile in_file_s (in_file);

  if (!in_file)
    {
      error ("audiowmark: error opening input playlist %s\n", in_name.c_str());
      return 1;
    }

  string out_name = out_dir + "/" + filename;
  FILE *out_file = fopen (out_name.c_str(), "w");
  ScopedFile out_file_s (out_file);

  if (!out_file)
    {
      error ("audiowmark: error opening output playlist %s\n", out_name.c_str());
      return 1;
    }

  struct Segment
  {
    string              name;
    map<string, string> vars;
  };
  vector<Segment> segments;
  char buffer[1024];
  int line = 1;
  const regex blank_re (R"(\s*(#.*)?)");
  while (fgets (buffer, 1024, in_file))
    {
      /* kill newline chars at end */
      int last = strlen (buffer) - 1;
      while (last > 0 && (buffer[last] == '\n' || buffer[last] == '\r'))
        buffer[last--] = 0;

      string s = buffer;

      std::smatch match;
      if (regex_match (s, blank_re))
        {
          /* blank line or comment */
          fprintf (out_file, "%s\n", s.c_str());
        }
      else
        {
          fprintf (out_file, "%s\n", s.c_str());
          Segment segment;
          segment.name = s;
          segments.push_back (segment);
        }
      line++;
    }
  size_t start_pos = 0;
  for (auto& segment : segments)
    {
      WavData out;
      Error err = ff_decode (in_dir + "/" + segment.name, out);
      if (err)
        {
          error ("audiowmark: hls: ff_decode failed: %s\n", err.message());
          return 1;
        }
      printf ("%d %zd\n", out.sample_rate(), out.n_values() / out.n_channels());
      segment.vars["start_pos"] = string_printf ("%zd", start_pos);
      start_pos += out.n_values() / out.n_channels();
    }
  for (size_t i = 0; i < segments.size(); i++)
    {
      TSWriter writer;

      if (i > 0)
        writer.append_file ("prev.ts", in_dir + "/" + segments[i - 1].name);
      if (i + 1 < segments.size())
        writer.append_file ("next.ts", in_dir + "/" + segments[i + 1].name);
      writer.append_vars ("vars", segments[i].vars);
      writer.process (in_dir + "/" + segments[i].name, out_dir + "/" + segments[i].name);
    }
  return 0;
}

int
hls_mark (const string& infile, const string& outfile, const string& bits)
{
  TSReader reader;

  Error err = reader.load (infile);
  if (err)
    {
      error ("hls_mark: %s\n", err.message());
      return 1;
    }

  WavData wav_data;
  err = ff_decode (infile, wav_data);
  if (err)
    {
      error ("hls_mark: %s\n", err.message());
      return 1;
    }

  for (auto entry : reader.entries())
    printf ("%s %zd\n", entry.filename.c_str(), entry.data.size());

  map<string, string> vars = reader.parse_vars ("vars");
  for (auto kv : vars)
    printf ("|| %s=%s\n", kv.first.c_str(), kv.second.c_str());

  size_t start_pos = atoi (vars["start_pos"].c_str());
  err = ff_encode (wav_data, outfile, start_pos);
  if (err)
    {
      error ("hls_mark: %s\n", err.message());
      return 1;
    }

  return 0;
}

int
main (int argc, char **argv)
{
  if (argc == 5 && strcmp (argv[1], "hls-embed-context") == 0)
    {
      printf ("hls-embed-context: in_dir=%s out_dir=%s m3u8=%s\n", argv[2], argv[3], argv[4]);
      return hls_embed_context (argv[2], argv[3], argv[4]);
    }
  else if (argc == 5 && strcmp (argv[1], "hls-mark") == 0)
    {
      return hls_mark (argv[2], argv[3], argv[4]);
    }
  else
    {
      error ("testhls: error parsing command line arguments\n");
    }
}

