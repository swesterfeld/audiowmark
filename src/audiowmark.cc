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
#include <math.h>
#include <string>
#include <random>
#include <algorithm>
#include <memory>

#include "wavdata.hh"
#include "utils.hh"
#include "random.hh"
#include "wmcommon.hh"
#include "shortcode.hh"
#include "hls.hh"

#include <assert.h>

#include "config.h"

using std::string;
using std::vector;
using std::min;
using std::max;

void
print_usage()
{
  printf ("usage: audiowmark <command> [ <args>... ]\n");
  printf ("\n");
  printf ("Commands:\n");
  printf ("  * create a watermarked wav file with a message\n");
  printf ("    audiowmark add <input_wav> <watermarked_wav> <message_hex>\n");
  printf ("\n");
  printf ("  * retrieve message\n");
  printf ("    audiowmark get <watermarked_wav>\n");
  printf ("\n");
  printf ("  * compare watermark message with expected message\n");
  printf ("    audiowmark cmp <watermarked_wav> <message_hex>\n");
  printf ("\n");
  printf ("  * generate 128-bit watermarking key, to be used with --key option\n");
  printf ("    audiowmark gen-key <key_file>\n");
  printf ("\n");
  printf ("Global options:\n");
  printf ("  --strength <s>        set watermark strength              [%.6g]\n", Params::water_delta * 1000);
  printf ("  --linear              disable non-linear bit storage\n");
  printf ("  --short <bits>        enable short payload mode\n");
  printf ("  --key <file>          load watermarking key from file\n");
  printf ("  -q, --quiet           disable information messages\n");
  printf ("\n");
  printf ("  --input-format raw    use raw stream as input\n");
  printf ("  --output-format raw   use raw stream as output\n");
  printf ("  --format raw          use raw stream as input and output\n");
  printf ("\n");
  printf ("The options to set the raw stream parameters (such as --raw-rate\n");
  printf ("or --raw-channels) are documented in the README file.\n");
}

Format
parse_format (const string& str)
{
  if (str == "raw")
    return Format::RAW;
  if (str == "auto")
    return Format::AUTO;
  error ("audiowmark: unsupported format '%s'\n", str.c_str());
  exit (1);
}

RawFormat::Endian
parse_endian (const string& str)
{
  if (str == "little")
    return RawFormat::Endian::LITTLE;
  if (str == "big")
    return RawFormat::Endian::BIG;
  error ("audiowmark: unsupported endianness '%s'\n", str.c_str());
  exit (1);
}

RawFormat::Encoding
parse_encoding (const string& str)
{
  if (str == "signed")
    return RawFormat::Encoding::SIGNED;
  if (str == "unsigned")
    return RawFormat::Encoding::UNSIGNED;
  error ("audiowmark: unsupported encoding '%s'\n", str.c_str());
  exit (1);
}

int
gentest (const string& infile, const string& outfile)
{
  printf ("generating test sample from '%s' to '%s'\n", infile.c_str(), outfile.c_str());

  WavData wav_data;
  Error err = wav_data.load (infile);
  if (err)
    {
      error ("audiowmark: error loading %s: %s\n", infile.c_str(), err.message());
      return 1;
    }
  const vector<float>& in_signal = wav_data.samples();
  vector<float> out_signal;

  /* 2:45 of audio - this is approximately the minimal amount of audio data required
   * for storing three separate watermarks with a 128-bit encoded message */
  const size_t offset = 0 * wav_data.n_channels() * wav_data.sample_rate();
  const size_t n_samples = 165 * wav_data.n_channels() * wav_data.sample_rate();
  if (in_signal.size() < (offset + n_samples))
    {
      error ("audiowmark: input file %s too short\n", infile.c_str());
      return 1;
    }
  for (size_t i = 0; i < n_samples; i++)
    {
      out_signal.push_back (in_signal[i + offset]);
    }
  WavData out_wav_data (out_signal, wav_data.n_channels(), wav_data.sample_rate(), wav_data.bit_depth());
  err = out_wav_data.save (outfile);
  if (err)
    {
      error ("audiowmark: error saving %s: %s\n", outfile.c_str(), err.message());
      return 1;
    }
  return 0;
}

int
cut_start (const string& infile, const string& outfile, const string& start_str)
{
  WavData wav_data;
  Error err = wav_data.load (infile);
  if (err)
    {
      error ("audiowmark: error loading %s: %s\n", infile.c_str(), err.message());
      return 1;
    }

  size_t start = atoi (start_str.c_str());

  const vector<float>& in_signal = wav_data.samples();
  vector<float> out_signal;
  for (size_t i = start * wav_data.n_channels(); i < in_signal.size(); i++)
    out_signal.push_back (in_signal[i]);

  WavData out_wav_data (out_signal, wav_data.n_channels(), wav_data.sample_rate(), wav_data.bit_depth());
  err = out_wav_data.save (outfile);
  if (err)
    {
      error ("audiowmark: error saving %s: %s\n", outfile.c_str(), err.message());
      return 1;
    }
  return 0;
}

int
test_subtract (const string& infile1, const string& infile2, const string& outfile)
{
  WavData in1_data;
  Error err = in1_data.load (infile1);
  if (err)
    {
      error ("audiowmark: error loading %s: %s\n", infile1.c_str(), err.message());
      return 1;
    }
  WavData in2_data;
  err = in2_data.load (infile2);
  if (err)
    {
      error ("audiowmark: error loading %s: %s\n", infile2.c_str(), err.message());
      return 1;
    }
  if (in1_data.n_values() != in2_data.n_values())
    {
      int64_t l1 = in1_data.n_values();
      int64_t l2 = in2_data.n_values();
      size_t  delta = std::abs (l1 - l2);
      warning ("audiowmark: size mismatch: %zd frames\n", delta / in1_data.n_channels());
      warning (" - %s frames: %zd\n", infile1.c_str(), in1_data.n_values() / in1_data.n_channels());
      warning (" - %s frames: %zd\n", infile2.c_str(), in2_data.n_values() / in2_data.n_channels());
    }
  assert (in1_data.n_channels() == in2_data.n_channels());

  const auto& in1_signal = in1_data.samples();
  const auto& in2_signal = in2_data.samples();
  size_t len = std::min (in1_data.n_values(), in2_data.n_values());
  vector<float> out_signal;
  for (size_t i = 0; i < len; i++)
    out_signal.push_back (in1_signal[i] - in2_signal[i]);

  WavData out_wav_data (out_signal, in1_data.n_channels(), in1_data.sample_rate(), in1_data.bit_depth());
  err = out_wav_data.save (outfile);
  if (err)
    {
      error ("audiowmark: error saving %s: %s\n", outfile.c_str(), err.message());
      return 1;
    }
  return 0;
}

int
test_snr (const string& orig_file, const string& wm_file)
{
  WavData orig_data;
  Error err = orig_data.load (orig_file);
  if (err)
    {
      error ("audiowmark: error loading %s: %s\n", orig_file.c_str(), err.message());
      return 1;
    }
  WavData wm_data;
  err = wm_data.load (wm_file);
  if (err)
    {
      error ("audiowmark: error loading %s: %s\n", wm_file.c_str(), err.message());
      return 1;
    }
  assert (orig_data.n_values() == wm_data.n_values());
  assert (orig_data.n_channels() == orig_data.n_channels());

  const auto& orig_signal = orig_data.samples();
  const auto& wm_signal = wm_data.samples();

  double snr_delta_power = 0;
  double snr_signal_power = 0;

  for (size_t i = 0; i < orig_signal.size(); i++)
    {
      const double orig  = orig_signal[i];                // original sample
      const double delta = orig_signal[i] - wm_signal[i]; // watermark

      snr_delta_power += delta * delta;
      snr_signal_power += orig * orig;
    }
  printf ("%f\n", 10 * log10 (snr_signal_power / snr_delta_power));
  return 0;
}

int
test_clip (const string& in_file, const string& out_file, int seed, int time_seconds)
{
  WavData in_data;
  Error err = in_data.load (in_file);
  if (err)
    {
      error ("audiowmark: error loading %s: %s\n", in_file.c_str(), err.message());
      return 1;
    }
  bool done = false;
  Random rng (seed, /* there is no stream for this test */ Random::Stream::data_up_down);
  size_t start_point, end_point;
  do
    {
      // this is unbiased only if 2 * block_size + time_seconds is smaller than overall file length
      const size_t values_per_block = (mark_sync_frame_count() + mark_data_frame_count()) * Params::frame_size * in_data.n_channels();
      start_point = 2 * values_per_block * (double(rng()) / UINT64_MAX);
      start_point /= in_data.n_channels();

      end_point = start_point + time_seconds * in_data.sample_rate();
      if (end_point < in_data.n_values() / in_data.n_channels())
        done = true;
    }
  while (!done);
  //printf ("%.3f %.3f\n", start_point / double (in_data.sample_rate()), end_point / double (in_data.sample_rate()));

  vector<float> out_signal (in_data.samples().begin() + start_point * in_data.n_channels(),
                            in_data.samples().begin() + end_point * in_data.n_channels());
  WavData out_wav_data (out_signal, in_data.n_channels(), in_data.sample_rate(), in_data.bit_depth());
  err = out_wav_data.save (out_file);
  if (err)
    {
      error ("audiowmark: error saving %s: %s\n", out_file.c_str(), err.message());
      return 1;
    }
  return 0;
}

int
gen_key (const string& outfile)
{
  FILE *f = fopen (outfile.c_str(), "w");
  if (!f)
    {
      error ("audiowmark: error writing to file %s\n", outfile.c_str());
      return 1;
    }
  fprintf (f, "# watermarking key for audiowmark\n\nkey %s\n", Random::gen_key().c_str());
  fclose (f);
  return 0;
}

class ArgParser
{
  vector<string> m_args;
  bool
  starts_with (const string& s, const string& start)
  {
    return s.substr (0, start.size()) == start;
  }
public:
  ArgParser (int argc, char **argv)
  {
    for (int i = 1; i < argc; i++)
      m_args.push_back (argv[i]);
  }
  bool
  parse_cmd (const string& cmd)
  {
    for (auto it = m_args.begin(); it != m_args.end(); it++)
      {
        if (!it->empty() && (*it)[0] != '-')
          {
            if (*it == cmd)
              {
                m_args.erase (it);
                return true;
              }
            else /* first positional arg is not cmd */
              {
                return false;
              }
          }
      }
    return false;
  }
  bool
  parse_opt (const string& option, string& out_s)
  {
    bool found_option = false;
    auto it = m_args.begin();
    while (it != m_args.end())
      {
        auto next_it = it + 1;
        if (*it == option && next_it != m_args.end())   /* --option foo */
          {
            out_s = *next_it;
            next_it = m_args.erase (it, it + 2);
            found_option = true;
          }
        else if (starts_with (*it, (option + "=")))   /* --option=foo */
          {
            out_s = it->substr (option.size() + 1);
            next_it = m_args.erase (it);
            found_option = true;
          }
        it = next_it;
      }
    return found_option;
  }
  bool
  parse_opt (const string& option, int& out_i)
  {
    string out_s;
    if (parse_opt (option, out_s))
      {
        out_i = atoi (out_s.c_str());
        return true;
      }
    return false;
  }
  bool
  parse_opt (const string& option, float& out_f)
  {
    string out_s;
    if (parse_opt (option, out_s))
      {
        out_f = atof (out_s.c_str());
        return true;
      }
    return false;
  }
  bool
  parse_opt (const string& option)
  {
    for (auto it = m_args.begin(); it != m_args.end(); it++)
      {
        if (*it == option) /* --option */
          {
            m_args.erase (it);
            return true;
          }
      }
    return false;
  }
  bool
  parse_args (size_t expected_count, vector<string>& out_args)
  {
    if (m_args.size() == expected_count)
      {
        out_args = m_args;
        return true;
      }
    return false;
  }
};

void
parse_shared_options (ArgParser& ap)
{
  int i;
  float f;
  string s;
  if (ap.parse_opt ("--strength", f))
    {
      Params::water_delta = f / 1000;
    }
  if (ap.parse_opt  ("--key", s))
    {
      Params::have_key++;
      Random::load_global_key (s);
    }
  if (ap.parse_opt ("--test-key", i))
    {
      Params::have_key++;
      Random::set_global_test_key (i);
    }
  if (ap.parse_opt ("--short", i))
    {
      Params::payload_size = i;
      if (!short_code_init (Params::payload_size))
        {
          error ("audiowmark: unsupported short payload size %zd\n", Params::payload_size);
          exit (1);
        }
      Params::payload_short = true;
    }
  ap.parse_opt ("--frames-per-bit", Params::frames_per_bit);
  if (ap.parse_opt ("--linear"))
    {
      Params::mix = false;
    }
  if (ap.parse_opt ("--quiet") || ap.parse_opt ("-q"))
    {
      set_log_level (Log::WARNING);
    }
  if (Params::have_key > 1)
    {
      error ("audiowmark: watermark key can at most be set once (--key / --test-key option)\n");
      exit (1);
    }
}

void
parse_add_options (ArgParser& ap)
{
  string s;
  int i;

  ap.parse_opt ("--set-input-label", Params::input_label);
  ap.parse_opt ("--set-output-label", Params::output_label);
  if (ap.parse_opt ("--snr"))
    {
      Params::snr = true;
    }
  if (ap.parse_opt ("--input-format", s))
    {
      Params::input_format = parse_format (s);
    }
  if (ap.parse_opt ("--output-format", s))
    {
      Params::output_format = parse_format (s);
    }
  if (ap.parse_opt ("--format", s))
    {
      Params::input_format = Params::output_format = parse_format (s);
    }
  if (ap.parse_opt ("--raw-input-bits", i))
    {
      Params::raw_input_format.set_bit_depth (i);
    }
  if (ap.parse_opt ("--raw-output-bits", i))
    {
      Params::raw_output_format.set_bit_depth (i);
    }
  if (ap.parse_opt ("--raw-bits", i))
    {
      Params::raw_input_format.set_bit_depth (i);
      Params::raw_output_format.set_bit_depth (i);
    }
  if (ap.parse_opt ( "--raw-input-endian", s))
    {
      auto e = parse_endian (s);
      Params::raw_input_format.set_endian (e);
    }
  if (ap.parse_opt ("--raw-output-endian", s))
    {
      auto e = parse_endian (s);
      Params::raw_output_format.set_endian (e);
    }
  if (ap.parse_opt ("--raw-endian", s))
    {
      auto e = parse_endian (s);
      Params::raw_input_format.set_endian (e);
      Params::raw_output_format.set_endian (e);
    }
  if (ap.parse_opt ("--raw-input-encoding", s))
    {
      auto e = parse_encoding (s);
      Params::raw_input_format.set_encoding (e);
    }
  if (ap.parse_opt ("--raw-output-encoding", s))
    {
      auto e = parse_encoding (s);
      Params::raw_output_format.set_encoding (e);
    }
  if (ap.parse_opt ("--raw-encoding", s))
    {
      auto e = parse_encoding (s);
      Params::raw_input_format.set_encoding (e);
      Params::raw_output_format.set_encoding (e);
    }
  if (ap.parse_opt ("--raw-channels", i))
    {
      Params::raw_input_format.set_channels (i);
      Params::raw_output_format.set_channels (i);
    }
  if (ap.parse_opt ("--raw-rate", i))
    {
      Params::raw_input_format.set_sample_rate (i);
      Params::raw_output_format.set_sample_rate (i);
    }
  if (ap.parse_opt ("--test-no-limiter"))
    {
      Params::test_no_limiter = true;
    }
}

void
parse_get_options (ArgParser& ap)
{
  ap.parse_opt ("--test-cut", Params::test_cut);
  ap.parse_opt ("--test-truncate", Params::test_truncate);

  if (ap.parse_opt ("--hard"))
    {
      Params::hard = true;
    }
  if (ap.parse_opt ("--test-no-sync"))
    {
      Params::test_no_sync = true;
    }
}

int
main (int argc, char **argv)
{
  ArgParser ap (argc, argv);
  vector<string> args;

  if (ap.parse_opt ("--help") || ap.parse_opt ("-h"))
    {
      print_usage();
      return 0;
    }
  if (ap.parse_opt ("--version") || ap.parse_opt ("-v"))
    {
      printf ("audiowmark %s\n", VERSION);
      return 0;
    }
  if (ap.parse_cmd ("hls-add"))
    {
      parse_shared_options (ap);

      ap.parse_opt ("--bit-rate", Params::hls_bit_rate);

      if (ap.parse_args (3, args))
        return hls_add (args[0], args[1], args[2]);
    }
  else if (ap.parse_cmd ("add"))
    {
      parse_shared_options (ap);
      parse_add_options (ap);

      if (ap.parse_args (3, args))
        return add_watermark (args[0], args[1], args[2]);
    }
  else if (ap.parse_cmd ("get"))
    {
      parse_shared_options (ap);
      parse_get_options (ap);

      if (ap.parse_args (1, args))
        return get_watermark (args[0], /* no ber */ "");
    }
  else if (ap.parse_cmd ("cmp"))
    {
      parse_shared_options (ap);
      parse_get_options (ap);

      if (ap.parse_args (2, args))
        return get_watermark (args[0], args[1]);
    }
  else if (ap.parse_cmd ("gen-key"))
    {
      if (ap.parse_args (1, args))
        return gen_key (args[0]);
    }
  else if (ap.parse_cmd ("gentest"))
    {
      if (ap.parse_args (2, args))
        return gentest (args[0], args[1]);
    }
  else if (ap.parse_cmd ("cut-start"))
    {
      if (ap.parse_args (3, args))
        return cut_start (args[0], args[1], args[2]);
    }
  else if (ap.parse_cmd ("test-subtract"))
    {
      if (ap.parse_args (3, args))
        return test_subtract (args[0], args[1], args[2]);
    }
  else if (ap.parse_cmd ("test-snr"))
    {
      if (ap.parse_args (2, args))
        return test_snr (args[0], args[1]);
    }
  else if (ap.parse_cmd ("test-clip"))
    {
      if (ap.parse_args (4, args))
        test_clip (args[0], args[1], atoi (args[2].c_str()), atoi (args[3].c_str()));
    }
  error ("audiowmark: error parsing commandline args (use audiowmark -h)\n");
  return 1;
}
