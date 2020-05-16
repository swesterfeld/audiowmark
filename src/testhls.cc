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
#include "wmcommon.hh"

using std::string;
using std::regex;
using std::vector;
using std::map;
using std::min;

Error
xsystem (const string& cmd)
{
  info ("+++ %s\n", cmd.c_str());
  int rc = system (cmd.c_str());
  int exit_status = WEXITSTATUS (rc);
  if (exit_status != 0)
    {
      error ("audiowmark: failed to execute command:\n%s\n", cmd.c_str());
      return Error (string_printf ("system failed / exit status %d", exit_status));
    }
  return Error::Code::NONE;
}

Error
ff_decode (const string& filename, const TSReader& reader, WavData& out_wav_data)
{
  FILE *tmp_file = tmpfile();
  ScopedFile tmp_file_s (tmp_file);
  string tmp_file_name = string_printf ("/dev/fd/%d", fileno (tmp_file));

  FILE *input_tmp_file = tmpfile();
  ScopedFile input_tmp_file_s (input_tmp_file);
  string input_tmp_file_name = string_printf ("/dev/fd/%d", fileno (input_tmp_file));

  /* build input file by concatenating previous ts and current ts */
  auto prev_ts = reader.find ("prev.ts");
  if (prev_ts)
    {
      // write previous ts
      size_t r = fwrite (prev_ts->data.data(), 1, prev_ts->data.size(), input_tmp_file);
      if (r != prev_ts->data.size())
        return Error (string_printf ("unable to write ff_decode:prev.ts to %s\n", input_tmp_file_name.c_str()));
    }
  // write current ts
  FILE *main = fopen (filename.c_str(), "r");
  ScopedFile main_s (main);
  int c;
  while ((c = fgetc (main)) >= 0)
    fputc (c, input_tmp_file);

  auto next_ts = reader.find ("next.ts");
  if (next_ts)
    {
      // write next ts
      size_t r = fwrite (next_ts->data.data(), 1, next_ts->data.size(), input_tmp_file);
      if (r != next_ts->data.size())
        return Error (string_printf ("unable to write ff_decode:next.ts to %s\n", input_tmp_file_name.c_str()));
    }
  fflush (input_tmp_file);
  string cmd = string_printf ("ffmpeg -v error -y -f mpegts -i %s -f wav %s", input_tmp_file_name.c_str(), tmp_file_name.c_str());
  Error err = xsystem (cmd.c_str());
  if (err)
    return err;

  err = out_wav_data.load (tmp_file_name);
  return err;
}

Error
ff_encode (const WavData& wav_data, const string& filename, size_t start_pos, size_t cut_start, size_t cut_end, double pts_start)
{
  FILE *tmp_file = tmpfile();
  ScopedFile tmp_file_s (tmp_file);
  string tmp_file_name = string_printf ("/dev/fd/%d", fileno (tmp_file));

  Error err = wav_data.save (tmp_file_name);

  string cmd = string_printf ("ffmpeg -v error -y -i %s -f mpegts -c:a aac '%s'", tmp_file_name.c_str(), filename.c_str());
  err = xsystem (cmd);
  if (err)
    return err;

  double length_s = double (wav_data.n_values()) / wav_data.n_channels() / wav_data.sample_rate();

  // cut_start_s is corrected down to avoid cutting one frame more or than intended
  double cut_start_s = cut_start / double (wav_data.sample_rate());
  double cut_end_s = cut_end / double (wav_data.sample_rate());
  cmd = string_printf ("ffmpeg -v error -y -i '%s' -ss %.6f -t %.6f -f mpegts -output_ts_offset %f -muxdelay 0 -muxpreload 0 -c copy '%s-tcpy'",
                       filename.c_str(), cut_start_s - 1. / wav_data.sample_rate(), length_s - (cut_start_s + cut_end_s), pts_start, filename.c_str());
  err = xsystem (cmd);
  if (err)
    return err;

  cmd = string_printf ("mv '%s-tcpy' '%s'", filename.c_str(), filename.c_str());
  xsystem (cmd);
  if (err)
    return err;

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
    size_t              size;
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
      Error err = ff_decode (in_dir + "/" + segment.name, /* FIXME: no context */ TSReader(), out);
      if (err)
        {
          error ("audiowmark: hls: ff_decode failed: %s\n", err.message());
          return 1;
        }
      printf ("%d %zd\n", out.sample_rate(), out.n_values() / out.n_channels());
      segment.size = out.n_values() / out.n_channels();

      /* obtain pts for first frame */
      string cmd = string_printf ("ffprobe -v 0 -show_entries packet=pts_time %s/%s -of compact=p=0:nk=1 | grep '^[0-9]'", in_dir.c_str(), segment.name.c_str());
      FILE *pts = popen (cmd.c_str(), "r");
      char buffer[1024];
      if (fgets (buffer, 1024, pts))
        {
          if (strlen (buffer) && buffer[strlen (buffer) - 1] == '\n')
            buffer[strlen (buffer) - 1] = 0;
          segment.vars["pts_start"] = buffer;
        }
      fclose (pts);

      segment.vars["start_pos"] = string_printf ("%zd", start_pos);
      segment.vars["size"] = string_printf ("%zd", segment.size);
      start_pos += segment.size;
    }

  /* fill out next/prev size fields */
  for (size_t i = 0; i < segments.size(); i++)
    {
      if (i > 0)
        segments[i].vars["prev_size"] = string_printf ("%zd", segments[i - 1].size);
      else
        segments[i].vars["prev_size"] = "0";

      if (i + 1 < segments.size())
        segments[i].vars["next_size"] = string_printf ("%zd", segments[i + 1].size);
      else
        segments[i].vars["next_size"] = "0";
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

class WDInputStream : public AudioInputStream
{
  WavData *wav_data;
  size_t   read_pos = 0;
public:
  WDInputStream (WavData *wav_data) :
    wav_data (wav_data)
  {
  }
  int
  bit_depth() const override
  {
    return wav_data->bit_depth();
  }
  int
  sample_rate() const override
  {
    return wav_data->sample_rate();
  }
  int
  n_channels() const override
  {
    return wav_data->n_channels();
  }
  size_t
  n_frames() const override
  {
    return wav_data->n_values() / wav_data->n_channels();
  }
  Error
  read_frames (std::vector<float>& samples, size_t count) override
  {
    size_t read_count = min (n_frames() - read_pos, count);

    const auto& wsamples = wav_data->samples();
    samples.assign (wsamples.begin() + read_pos * n_channels(), wsamples.begin() + (read_pos + read_count) * n_channels());

    read_pos += read_count;

    return Error::Code::NONE;
  }
};

class WDOutputStream : public AudioOutputStream
{
  WavData *wav_data;
  vector<float> samples;
public:
  WDOutputStream (WavData *wav_data) :
    wav_data (wav_data)
  {
  }
  int
  bit_depth() const override
  {
    return wav_data->bit_depth();
  }
  int
  sample_rate() const override
  {
    return wav_data->sample_rate();
  }
  int
  n_channels() const override
  {
    return wav_data->n_channels();
  }
  Error
  write_frames (const std::vector<float>& frames)
  {
    samples.insert (samples.end(), frames.begin(), frames.end());
    return Error::Code::NONE;
  }
  Error
  close()
  {
    wav_data->set_samples (samples); // only do this once at end for performance reasons
    return Error::Code::NONE;
  }
};

int
mark_zexpand (WavData& wav_data, size_t zero_frames, const string& bits)
{
  WDInputStream in_stream (&wav_data);

  WavData wav_data_out ({ /* no samples */ }, wav_data.n_channels(), wav_data.sample_rate(), wav_data.bit_depth());
  WDOutputStream out_stream (&wav_data_out);

  int rc = add_stream_watermark (&in_stream, &out_stream, bits, zero_frames);
  if (rc != 0)
    return rc;

  wav_data.set_samples (wav_data_out.samples());

  return 0;
}

int
hls_mark (const string& infile, const string& outfile, const string& bits)
{
  double start_time = get_time();

  TSReader reader;

  Error err = reader.load (infile);
  if (err)
    {
      error ("hls_mark: %s\n", err.message());
      return 1;
    }

  WavData wav_data;
  err = ff_decode (infile, reader, wav_data);
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
  size_t prev_size = atoi (vars["prev_size"].c_str());
  size_t next_size = atoi (vars["next_size"].c_str());
  double pts_start = atof (vars["pts_start"].c_str());
  size_t next_ctx = min<size_t> (1024 * 3, next_size);
  size_t prev_ctx = min<size_t> (1024 * 3, prev_size);

  int zrc = mark_zexpand (wav_data, start_pos - prev_size, bits);
  if (zrc != 0)
    return zrc;

  /* erase extra samples caused by concatting with prev.ts */
  auto samples = wav_data.samples();
  samples.erase (samples.begin(), samples.begin() + (prev_size - prev_ctx) * wav_data.n_channels());
  samples.erase (samples.end() - (next_size - next_ctx) * wav_data.n_channels(), samples.end());
  wav_data.set_samples (samples);

  err = ff_encode (wav_data, outfile, start_pos, start_pos == 0 ? 1024 : prev_ctx, next_ctx, pts_start);
  if (err)
    {
      error ("hls_mark: %s\n", err.message());
      return 1;
    }

  double end_time = get_time();
  printf ("hls_time %f %f\n", start_pos / double (wav_data.sample_rate()), (end_time - start_time) * 1000 /* ms */);

  return 0;
}

int
test_seek (const string& in, const string& out, int pos, const string& bits)
{
  vector<float> samples;
  WavData wav_data;
  Error err = wav_data.load (in);
  if (err)
    {
      error ("load error: %s\n", err.message());
      return 1;
    }

  samples = wav_data.samples();
  samples.erase (samples.begin(), samples.begin() + pos * wav_data.n_channels());
  wav_data.set_samples (samples);

  int rc = mark_zexpand (wav_data, pos, bits);
  if (rc != 0)
    {
      return rc;
    }

  samples = wav_data.samples();
  samples.insert (samples.begin(), pos * wav_data.n_channels(), 0);
  wav_data.set_samples (samples);

  err = wav_data.save (out);
  if (err)
    {
      error ("save error: %s\n", err.message());
      return 1;
    }
  return 0;
}

int
seek_perf (int sample_rate, int seconds)
{
  vector<float> samples (100);
  WavData wav_data (samples, 2, sample_rate, 16);

  double start_time = get_time();

  int rc = mark_zexpand (wav_data, seconds * sample_rate, "0c");
  if (rc != 0)
    return rc;

  double end_time = get_time();

  printf ("\n\n");
  printf ("total time %7.3f sec\n", end_time - start_time);
  printf ("per second %7.3f ms\n", (end_time - start_time) / seconds * 1000);

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
  else if (argc == 6 && strcmp (argv[1], "test-seek") == 0)
    {
      return test_seek (argv[2], argv[3], atoi (argv[4]), argv[5]);
    }
  else if (argc == 4 && strcmp (argv[1], "seek-perf") == 0)
    {
      return seek_perf (atoi (argv[2]), atoi (argv[3]));
    }
  else
    {
      error ("testhls: error parsing command line arguments\n");
      return 1;
    }
}

