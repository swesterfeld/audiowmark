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

#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <unistd.h>

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

string
args2string (const vector<string>& args)
{
  string result;
  bool first = true;

  for (auto a : args)
    {
      if (!first)
        result += " ";
      first = false;
      result += a;
    }
  return result;
}

Error
run (const vector<string>& args, vector<string> *pipe_out = nullptr)
{
  auto report_error = [=] { error ("audiowmark: failed to execute %s\n", args2string (args).c_str()); };

  char *argv[args.size() + 1];
  for (size_t i = 0; i < args.size(); i++)
    argv[i] = (char *) args[i].c_str();
  argv[args.size()] = nullptr;

  int pipe_fds[2];
  if (pipe_out)
    {
      if (pipe (pipe_fds) == -1)
        {
          report_error();
          return Error ("pipe() failed");
        }
    }
  pid_t pid = fork();
  if (pid < 0)
    {
      if (pipe_out)
        {
          close (pipe_fds[0]);
          close (pipe_fds[1]);
        }
      report_error();
      return Error ("fork() failed");
    }
  if (pid == 0) /* child process */
    {
      if (pipe_out)
        {
          // replace stdout with pipe
          if (dup2 (pipe_fds[1], STDOUT_FILENO) == -1)
            {
              perror ("audiowmark: dup2() failed");
              exit (127);
            }

          // close remaining pipe fds
          close (pipe_fds[0]);
          close (pipe_fds[1]);
        }
      execvp (argv[0], argv);
      perror ("audiowmark: execvp() failed");

      // should not be reached in normal operation, so exec failed
      exit (127);
    }

  /* parent process */
  if (pipe_out) /* capture child stdout */
    {
      close (pipe_fds[1]); // close pipe write fd

      FILE *f = fdopen (pipe_fds[0], "r");
      if (!f)
        {
          close (pipe_fds[0]);
          report_error();
          return Error ("fdopen() pipe failed");
        }
      char buffer[1024];
      while (fgets (buffer, 1024, f))
        {
          if (strlen (buffer) && buffer[strlen (buffer) - 1] == '\n')
            buffer[strlen (buffer) - 1] = 0;
          if (pipe_out)
            pipe_out->push_back (buffer);
        }
      fclose (f);        // close pipe read fd
    }
  int status;
  pid_t exited = waitpid (pid, &status, 0);
  if (exited < 0)
    {
      report_error();
      return Error ("waitpid() failed");
    }
  if (WIFEXITED (status))
    {
      int exit_status = WEXITSTATUS (status);
      if (exit_status != 0)
        {
          report_error();
          return Error (string_printf ("subprocess failed / exit status %d", exit_status));
        }
    }
  else
    {
      report_error();
      return Error ("child didn't exit normally");
    }
  return Error::Code::NONE;
}

Error
ff_decode (const string& filename, WavData& out_wav_data)
{
  FILE *tmp_file = tmpfile();
  ScopedFile tmp_file_s (tmp_file);
  string tmp_file_name = string_printf ("/dev/fd/%d", fileno (tmp_file));

  if (!tmp_file)
    return Error ("failed to create temp file");

  Error err = run ({"ffmpeg", "-v", "error", "-y", "-f",  "mpegts", "-i", filename, "-f", "wav", tmp_file_name});
  if (err)
    return err;

  err = out_wav_data.load (tmp_file_name);
  return err;
}

int
hls_add (const string& infile, const string& outfile, const string& bits)
{
  TSReader reader;

  Error err = reader.load (infile);
  if (err)
    {
      error ("hls_mark: %s\n", err.message());
      return 1;
    }

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

  int wm_rc = add_stream_watermark (&in_stream, &out_stream, bits, start_pos - prev_size);
  if (wm_rc != 0)
    return wm_rc;

  info ("AAC Bitrate:  %d\n", bit_rate);
  return 0;
}

Error
bit_rate_from_m3u8 (const string& m3u8, const WavData& wav_data, int& bit_rate)
{
  FILE *tmp_file = tmpfile();
  ScopedFile tmp_file_s (tmp_file);
  string tmp_file_name = string_printf ("/dev/fd/%d", fileno (tmp_file));

  if (!tmp_file)
    return Error ("failed to create temp file");

  Error err = run ({"ffmpeg", "-v", "error", "-y", "-i", m3u8, "-c:a", "copy", "-f", "adts", tmp_file_name});
  if (err)
    return err;

  struct stat stat_buf;
  if (stat (tmp_file_name.c_str(), &stat_buf) != 0)
    {
      return Error (string_printf ("failed to stat temporary aac file: %s", strerror (errno)));
    }
  double seconds = double (wav_data.n_frames()) / wav_data.sample_rate();
  bit_rate = stat_buf.st_size / seconds * 8;
  return Error::Code::NONE;
}

Error
validate_input_segment (const string& filename)
{
  vector<string> format_out;
  Error err = run ({"ffprobe", "-v", "error", "-print_format", "compact", "-show_streams", filename}, &format_out);
  if (err)
    {
      error ("audiowmark: hls: failed to validate input file: %s\n", filename.c_str());
      return err;
    }
  for (auto o : format_out)
    {
      /* parse assignments stream|index=0|codec_name=aac|... */
      map<string, string> params;

      string key, value;
      bool in_key = true;
      for (char c : '|' + o + '|')
        {
          if (c == '=')
            {
              in_key = false;
            }
          else if (c == '|')
            {
              params[key] = value;
              in_key = true;
              key = "";
              value = "";
            }
          else
            {
              if (in_key)
                key += c;
              else
                value += c;
            }
        }

      /* now the actual validation */
      if (atoi (params["index"].c_str()) != 0)
        {
          error ("audiowmark: hls segment '%s' contains more than one stream\n", filename.c_str());
          return Error ("failed to validate input segment");
        }
      if (params["codec_name"] != "aac")
        {
          error ("audiowmark hls segment '%s' is not encoded using AAC\n", filename.c_str());
          return Error ("failed to validate input segment");
        }
    }
  return Error::Code::NONE;
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
  for (auto& segment : segments)
    {
      Error err = validate_input_segment (in_dir + "/" + segment.name);
      if (err)
        {
          error ("audiowmark: hls: %s\n", err.message());
          return 1;
        }
    }

  /* find bitrate for AAC encoder */
  int bit_rate = 0;
  if (!Params::hls_bit_rate)
    {
      err = bit_rate_from_m3u8 (in_dir + "/" + filename, audio_master_data, bit_rate);
      if (err)
        {
          error ("audiowmark: bit-rate detection failed: %s\n", err.message());
          return 1;
         }
      info ("AAC Bitrate:  %d (detected)\n", bit_rate);
    }
  else
    {
      bit_rate = Params::hls_bit_rate;
      info ("AAC Bitrate:  %d\n", bit_rate);
    }

  info ("Segments:     %zd\n", segments.size());
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
      segment.size = out.n_values() / out.n_channels();

      if ((segment.size % 1024) != 0)
        {
          error ("audiowmark: hls input segments need 1024-sample alignment (due to AAC)\n");
          return 1;
        }

      /* obtain pts for first frame */
      vector<string> pts_out;
      err = run ({"ffprobe", "-v", "0", "-show_entries", "packet=pts_time", in_dir + "/" + segment.name, "-of", "compact=p=0:nk=1"}, &pts_out);
      if (err)
        {
          error ("audiowmark: hls: initial pts detection failed: %s\n", err.message());
          return 1;
        }
      for (auto o : pts_out)
        {
          if (!o.empty())
            {
              segment.vars["pts_start"] = o;
              break;
            }
        }

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

      err = writer.process (in_dir + "/" + segment.name, out_dir + "/" + segment.name);
      if (err)
        {
          error ("audiowmark: processing hls segment %s failed: %s\n", segment.name.c_str(), err.message());
          return 1;
        }

      /* start position for the next segment */
      start_pos += segment.size;
    }
  int orig_seconds = start_pos / audio_master_data.sample_rate();
  info ("Time:         %d:%02d\n", orig_seconds / 60, orig_seconds % 60);
  return 0;
}


