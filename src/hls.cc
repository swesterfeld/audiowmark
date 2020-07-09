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

#include <string>
#include <regex>

#include "utils.hh"
#include "mpegts.hh"
#include "sfinputstream.hh"
#include "hlsoutputstream.hh"
#include "sfoutputstream.hh"
#include "wmcommon.hh"
#include "wavdata.hh"

using std::string;
using std::vector;
using std::regex;
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
ff_decode (const string& filename, WavData& out_wav_data)
{
  FILE *tmp_file = tmpfile();
  ScopedFile tmp_file_s (tmp_file);
  string tmp_file_name = string_printf ("/dev/fd/%d", fileno (tmp_file));

  FILE *input_tmp_file = tmpfile();
  ScopedFile input_tmp_file_s (input_tmp_file);
  string input_tmp_file_name = string_printf ("/dev/fd/%d", fileno (input_tmp_file));

  // write current ts
  FILE *main = fopen (filename.c_str(), "r");
  ScopedFile main_s (main);
  int c;
  while ((c = fgetc (main)) >= 0)
    fputc (c, input_tmp_file);

  fflush (input_tmp_file);
  string cmd = string_printf ("ffmpeg -v error -y -f mpegts -i %s -f wav %s", input_tmp_file_name.c_str(), tmp_file_name.c_str());
  Error err = xsystem (cmd.c_str());
  if (err)
    return err;

  err = out_wav_data.load (tmp_file_name);
  return err;
}

int
hls_add (const string& infile, const string& outfile, const string& bits)
{
  double start_time = get_time();

  TSReader reader;

  Error err = reader.load (infile);
  if (err)
    {
      error ("hls_mark: %s\n", err.message());
      return 1;
    }
  info ("hls_elapsed_load %f\n", (get_time() - start_time) * 1000 /* ms */);
  double start_time1 = get_time();

  const TSReader::Entry *full_flac = reader.find ("full.flac");
  if (!full_flac)
    {
      error ("hls_mark: no embedded context found in %s\n", infile.c_str());
      return 1;
    }

  SFInputStream in_stream;
  err = in_stream.open (&full_flac->data);
  if (err)
    {
      error ("hls_mark: %s\n", err.message());
      return 1;
    }

  for (auto entry : reader.entries())
    info ("%s %zd\n", entry.filename.c_str(), entry.data.size());

  map<string, string> vars = reader.parse_vars ("vars");
  for (auto kv : vars)
    info ("|| %s=%s\n", kv.first.c_str(), kv.second.c_str());

  size_t start_pos = atoi (vars["start_pos"].c_str());
  size_t prev_size = atoi (vars["prev_size"].c_str());
  size_t next_size = atoi (vars["next_size"].c_str());
  size_t size      = atoi (vars["size"].c_str());
  double pts_start = atof (vars["pts_start"].c_str());
  int    bit_rate = atoi (vars["bit_rate"].c_str());
  size_t prev_ctx = min<size_t> (1024 * 3, prev_size);

  info ("hls_time_elapsed_decode %f\n", (get_time() - start_time1) * 1000 /* ms */);
  start_time1 = get_time();

  if (Params::hls_bit_rate)  // command line option overrides vars bit-rate
    bit_rate = Params::hls_bit_rate;

  HLSOutputStream out_stream (in_stream.n_channels(), in_stream.sample_rate(), in_stream.bit_depth());

  out_stream.set_bit_rate (bit_rate);
  info ("n_frames = %zd\n", in_stream.n_frames() - prev_size - next_size);
  const size_t shift = 1024;
  const size_t cut_aac_frames = (prev_ctx + shift) / 1024;
  const size_t delete_input_start = prev_size - prev_ctx;
  const size_t keep_aac_frames = size / 1024;

  err = out_stream.open (outfile, cut_aac_frames, keep_aac_frames, pts_start, delete_input_start);
  if (err)
    {
      error ("audiowmark: error opening HLS output stream %s: %s\n", outfile.c_str(), err.message());
      return 1;
    }

  int zrc = add_stream_watermark (&in_stream, &out_stream, bits, start_pos - prev_size);
  if (zrc != 0)
    {
      info ("hls_time_abort_enc %f\n", (get_time() - start_time1) * 1000 /* ms */);

      double end_time = get_time();
      info ("hls_time_abort %f %f\n", start_pos / double (out_stream.sample_rate()), (end_time - start_time) * 1000 /* ms */);
      return zrc;
    }
  info ("AAC Bitrate:  %d\n", bit_rate);

  info ("hls_time_elapsed_aac_enc %f\n", (get_time() - start_time1) * 1000 /* ms */);

  double end_time = get_time();
  info ("hls_time %f %f\n", start_pos / double (out_stream.sample_rate()), (end_time - start_time) * 1000 /* ms */);

  return 0;
}

int
hls_prepare (const string& in_dir, const string& out_dir, const string& filename, const string& audio_master)
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

  WavData audio_master_data;
  Error err = audio_master_data.load (audio_master);
  if (err)
    {
      error ("audiowmark: failed to load audio master: %s\n", audio_master.c_str());
      return 1;
    }

  const int bit_rate = Params::hls_bit_rate;
  info ("AAC Bitrate:  %d\n", bit_rate);

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
      Error err = ff_decode (in_dir + "/" + segment.name, out);
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

      /* store 3 seconds of the context before this segment and after this segment (if available) */
      const size_t ctx_3sec = 3 * out.sample_rate();
      const size_t prev_size = min<size_t> (start_pos, ctx_3sec);
      const size_t next_size = min<size_t> (audio_master_data.n_frames() - (segment.size + start_pos), ctx_3sec);

      segment.vars["start_pos"] = string_printf ("%zd", start_pos);
      segment.vars["size"] = string_printf ("%zd", segment.size);
      segment.vars["prev_size"] = string_printf ("%zd", prev_size);
      segment.vars["next_size"] = string_printf ("%zd", next_size);
      segment.vars["bit_rate"] = string_printf ("%d", bit_rate);

      /* write audio segment with context */
      const size_t start_point = start_pos - prev_size;
      const size_t end_point = start_point + prev_size + segment.size + next_size;

      vector<float> out_signal (audio_master_data.samples().begin() + start_point * audio_master_data.n_channels(),
                                audio_master_data.samples().begin() + end_point * audio_master_data.n_channels());

      vector<unsigned char> full_flac_mem;
      SFOutputStream out_stream;
      err = out_stream.open (&full_flac_mem,
                             audio_master_data.n_channels(), audio_master_data.sample_rate(), audio_master_data.bit_depth(),
                             SFOutputStream::OutFormat::FLAC);
      if (err)
        {
          error ("audiowmark: hls: open context flac failed: %s\n", err.message());
          return 1;
        }

      err = out_stream.write_frames (out_signal);
      if (err)
        {
          error ("audiowmark: hls: write context flac failed: %s\n", err.message());
          return 1;
        }

      err = out_stream.close();
      if (err)
        {
          error ("audiowmark: hls: close context flac failed: %s\n", err.message());
          return 1;
        }

      /* store everything we need in a mpegts file */
      TSWriter writer;

      writer.append_data ("full.flac", full_flac_mem);
      writer.append_vars ("vars", segment.vars);
      writer.process (in_dir + "/" + segment.name, out_dir + "/" + segment.name);

      /* start position for the next segment */
      start_pos += segment.size;
    }
  return 0;
}


