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
#include "sfoutputstream.hh"
#include "wmcommon.hh"
#include "wavdata.hh"

#include "config.h"

using std::string;
using std::vector;
using std::regex;
using std::map;
using std::min;

#if !HAVE_FFMPEG
int
hls_prepare (const string& in_dir, const string& out_dir, const string& filename, const string& audio_master)
{
  error ("audiowmark: hls support is not available in this build of audiowmark\n");
  return 1;
}

int
hls_add (const string& infile, const string& outfile, const string& bits)
{
  error ("audiowmark: hls support is not available in this build of audiowmark\n");
  return 1;
}
#else

#include "hlsoutputstream.hh"

static bool
file_exists (const string& filename)
{
  struct stat st;

  if (stat (filename.c_str(), &st) == 0)
    {
      return S_ISREG (st.st_mode);
    }
  return false;
}

static string
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

static Error
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
      error ("hls: %s\n", err.message());
      return 1;
    }

  const TSReader::Entry *full_flac = reader.find ("full.flac");
  if (!full_flac)
    {
      error ("hls: no embedded context found in %s\n", infile.c_str());
      return 1;
    }

  SFInputStream in_stream;
  err = in_stream.open (&full_flac->data);
  if (err)
    {
      error ("hls: %s\n", err.message());
      return 1;
    }

  map<string, string> vars = reader.parse_vars ("vars");
  bool missing_vars = false;

  auto get_var = [&] (const std::string& var) {
    auto it = vars.find (var);
    if (it == vars.end())
      {
        error ("audiowmark: hls segment is missing value for required variable '%s'\n", var.c_str());
        missing_vars = true;
        return "";
      }
    else
      return it->second.c_str();
  };
  size_t start_pos = atoi (get_var ("start_pos"));
  size_t prev_size = atoi (get_var ("prev_size"));
  size_t size      = atoi (get_var ("size"));
  double pts_start = atof (get_var ("pts_start"));
  int    bit_rate  = atoi (get_var ("bit_rate"));
  size_t prev_ctx  = min<size_t> (1024 * 3, prev_size);

  string channel_layout = get_var ("channel_layout");

  if (missing_vars)
    return 1;

  if (Params::hls_bit_rate)  // command line option overrides vars bit-rate
    bit_rate = Params::hls_bit_rate;

  HLSOutputStream out_stream (in_stream.n_channels(), in_stream.sample_rate(), in_stream.bit_depth());

  out_stream.set_bit_rate (bit_rate);
  out_stream.set_channel_layout (channel_layout);

  /* ffmpeg aac encode adds one frame of latency - it would be possible to compensate for this
   * by setting shift = 1024, but it can also be done by adjusting the presentation timestamp
   */
  const size_t shift = 0;
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
load_audio_master (const string& filename, WavData& audio_master_data)
{
  FILE *tmp_file = tmpfile();
  ScopedFile tmp_file_s (tmp_file);
  string tmp_file_name = string_printf ("/dev/fd/%d", fileno (tmp_file));

  if (!tmp_file)
    return Error ("failed to create temp file");

  /* extract wav */
  Error err = run ({"ffmpeg", "-v", "error", "-y", "-i", filename, "-f", "wav", tmp_file_name});
  if (err)
    return err;

  err = audio_master_data.load (tmp_file_name);
  if (err)
    return err;

  return Error::Code::NONE;
}

Error
probe_input_segment (const string& filename, map<string, string>& params)
{
  TSReader reader;

  Error err = reader.load (filename);
  if (err)
    {
      error ("audiowmark: hls: failed to read mpegts input file: %s\n", filename.c_str());
      return err;
    }

  if (reader.entries().size())
    {
      error ("audiowmark: hls: file appears to be already prepared: %s\n", filename.c_str());
      return Error ("input for hls-prepare must not contain context");
    }

  vector<string> format_out;
  err = run ({"ffprobe", "-v", "error", "-print_format", "compact", "-show_streams", filename}, &format_out);
  if (err)
    {
      error ("audiowmark: hls: failed to validate input file: %s\n", filename.c_str());
      return err;
    }
  for (auto o : format_out)
    {
      /* parse assignments stream|index=0|codec_name=aac|... */
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

  int mkret = mkdir (out_dir.c_str(), 0755);
  if (mkret == -1 && errno != EEXIST)
    {
      error ("audiowmark: unable to create directory %s: %s\n", out_dir.c_str(), strerror (errno));
      return 1;
    }

  string out_name = out_dir + "/" + filename;
  if (file_exists (out_name))
    {
      error ("audiowmark: output file already exists: %s\n", out_name.c_str());
      return 1;
    }
  FILE *out_file = fopen (out_name.c_str(), "w");
  ScopedFile out_file_s (out_file);

  if (!out_file)
    {
      error ("audiowmark: error opening output playlist %s\n", out_name.c_str());
      return 1;
    }

  WavData audio_master_data;
  Error err = load_audio_master (audio_master, audio_master_data);
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
      map<string, string> params;
      string segname = in_dir + "/" + segment.name;

      Error err = probe_input_segment (segname, params);
      if (err)
        {
          error ("audiowmark: hls: %s\n", err.message());
          return 1;
        }
      /* validate input segment */
      if (atoi (params["index"].c_str()) != 0)
        {
          error ("audiowmark: hls segment '%s' contains more than one stream\n", segname.c_str());
          return 1;
        }
      if (params["codec_name"] != "aac")
        {
          error ("audiowmark: hls segment '%s' is not encoded using AAC\n", segname.c_str());
          return 1;
        }

      /* get segment parameters */
      if (params["channel_layout"].empty())
        {
          error ("audiowmark: hls segment '%s' has no channel_layout entry\n", segname.c_str());
          return 1;
        }
      segment.vars["channel_layout"] = params["channel_layout"];

      /* get start pts */
      if (params["start_time"].empty())
        {
          error ("audiowmark: hls segment '%s' has no start_time entry\n", segname.c_str());
          return 1;
        }
      segment.vars["pts_start"] = params["start_time"];
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

      /* store 3 seconds of the context before this segment and after this segment (if available) */
      const size_t ctx_3sec = 3 * out.sample_rate();
      const size_t prev_size = min<size_t> (start_pos, ctx_3sec);
      const size_t segment_size_with_ctx = prev_size + segment.size + ctx_3sec;

      segment.vars["start_pos"] = string_printf ("%zd", start_pos);
      segment.vars["size"] = string_printf ("%zd", segment.size);
      segment.vars["prev_size"] = string_printf ("%zd", prev_size);
      segment.vars["bit_rate"] = string_printf ("%d", bit_rate);

      /* write audio segment with context */
      const size_t start_point = min (start_pos - prev_size, audio_master_data.n_frames());
      const size_t end_point = min (start_point + segment_size_with_ctx, audio_master_data.n_frames());

      vector<float> out_signal (audio_master_data.samples().begin() + start_point * audio_master_data.n_channels(),
                                audio_master_data.samples().begin() + end_point * audio_master_data.n_channels());

      // append zeros if audio master is too short to provide segment with context
      out_signal.resize (segment_size_with_ctx * audio_master_data.n_channels());

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

      string out_segment = out_dir + "/" + segment.name;
      if (file_exists (out_segment))
        {
          error ("audiowmark: output file already exists: %s\n", out_segment.c_str());
          return 1;
        }
      err = writer.process (in_dir + "/" + segment.name, out_segment);
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
#endif
