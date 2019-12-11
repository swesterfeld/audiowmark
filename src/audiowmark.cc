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
  printf ("  --key <file>          load watermarking key from file\n");
  printf ("  -q, --quiet           disable information messages\n");
}

static bool
check_arg (uint         argc,
           char        *argv[],
           uint        *nth,
           const char  *opt,		    /* for example: --foo */
           const char **opt_arg = nullptr)  /* if foo needs an argument, pass a pointer to get the argument */
{
  assert (opt != nullptr);
  assert (*nth < argc);

  const char *arg = argv[*nth];
  if (!arg)
    return false;

  uint opt_len = strlen (opt);
  if (strcmp (arg, opt) == 0)
    {
      if (opt_arg && *nth + 1 < argc)	  /* match foo option with argument: --foo bar */
        {
          argv[(*nth)++] = nullptr;
          *opt_arg = argv[*nth];
          argv[*nth] = nullptr;
          return true;
        }
      else if (!opt_arg)		  /* match foo option without argument: --foo */
        {
          argv[*nth] = nullptr;
          return true;
        }
      /* fall through to error message */
    }
  else if (strncmp (arg, opt, opt_len) == 0 && arg[opt_len] == '=')
    {
      if (opt_arg)			  /* match foo option with argument: --foo=bar */
        {
          *opt_arg = arg + opt_len + 1;
          argv[*nth] = nullptr;
          return true;
        }
      /* fall through to error message */
    }
  else
    return false;

  print_usage();
  exit (1);
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

void
parse_options (int   *argc_p,
               char **argv_p[])
{
  uint argc = *argc_p;
  char **argv = *argv_p;
  unsigned int i, e;

  for (i = 1; i < argc; i++)
    {
      const char *opt_arg;
      if (strcmp (argv[i], "--help") == 0 ||
          strcmp (argv[i], "-h") == 0)
	{
	  print_usage();
	  exit (0);
	}
      else if (strcmp (argv[i], "--version") == 0 || strcmp (argv[i], "-v") == 0)
	{
	  printf ("audiowmark %s\n", VERSION);
	  exit (0);
	}
      else if (check_arg (argc, argv, &i, "--frames-per-bit", &opt_arg))
	{
          Params::frames_per_bit = atoi (opt_arg);
	}
      else if (check_arg (argc, argv, &i, "--strength", &opt_arg))
	{
          Params::water_delta = atof (opt_arg) / 1000;
	}
      else if (check_arg (argc, argv, &i, "--linear"))
	{
          Params::mix = false;
	}
      else if (check_arg (argc, argv, &i, "--hard"))
	{
          Params::hard = true;
	}
      else if (check_arg (argc, argv, &i, "--snr"))
        {
          Params::snr = true;
        }
      else if (check_arg (argc, argv, &i, "--test-key", &opt_arg))
	{
          Params::have_key++;
          Random::set_global_test_key (atoi (opt_arg));
	}
      else if (check_arg (argc, argv, &i, "--key", &opt_arg))
        {
          Params::have_key++;
          Random::load_global_key (opt_arg);
        }
      else if (check_arg (argc, argv, &i, "--test-cut", &opt_arg))
	{
          Params::test_cut = atoi (opt_arg);
	}
      else if (check_arg (argc, argv, &i, "--test-no-sync"))
        {
          Params::test_no_sync = true;
        }
      else if (check_arg (argc, argv, &i, "--test-truncate", &opt_arg))
	{
          Params::test_truncate = atoi (opt_arg);
	}
      else if (check_arg (argc, argv, &i, "--input-format", &opt_arg))
        {
          Params::input_format = parse_format (opt_arg);
        }
      else if (check_arg (argc, argv, &i, "--output-format", &opt_arg))
        {
          Params::output_format = parse_format (opt_arg);
        }
      else if (check_arg (argc, argv, &i, "--format", &opt_arg))
        {
          Params::input_format = Params::output_format = parse_format (opt_arg);
        }
      else if (check_arg (argc, argv, &i, "--raw-input-bits", &opt_arg))
        {
          int b = atoi (opt_arg);
          Params::raw_input_format.set_bit_depth (b);
        }
      else if (check_arg (argc, argv, &i, "--raw-output-bits", &opt_arg))
        {
          int b = atoi (opt_arg);
          Params::raw_output_format.set_bit_depth (b);
        }
      else if (check_arg (argc, argv, &i, "--raw-bits", &opt_arg))
        {
          int b = atoi (opt_arg);
          Params::raw_input_format.set_bit_depth (b);
          Params::raw_output_format.set_bit_depth (b);
        }
      else if (check_arg (argc, argv, &i, "--raw-input-endian", &opt_arg))
        {
          auto e = parse_endian (opt_arg);
          Params::raw_input_format.set_endian (e);
        }
      else if (check_arg (argc, argv, &i, "--raw-output-endian", &opt_arg))
        {
          auto e = parse_endian (opt_arg);
          Params::raw_output_format.set_endian (e);
        }
      else if (check_arg (argc, argv, &i, "--raw-endian", &opt_arg))
        {
          auto e = parse_endian (opt_arg);
          Params::raw_input_format.set_endian (e);
          Params::raw_output_format.set_endian (e);
        }
      else if (check_arg (argc, argv, &i, "--raw-input-encoding", &opt_arg))
        {
          auto e = parse_encoding (opt_arg);
          Params::raw_input_format.set_encoding (e);
        }
      else if (check_arg (argc, argv, &i, "--raw-output-encoding", &opt_arg))
        {
          auto e = parse_encoding (opt_arg);
          Params::raw_output_format.set_encoding (e);
        }
      else if (check_arg (argc, argv, &i, "--raw-encoding", &opt_arg))
        {
          auto e = parse_encoding (opt_arg);
          Params::raw_input_format.set_encoding (e);
          Params::raw_output_format.set_encoding (e);
        }
      else if (check_arg (argc, argv, &i, "--raw-channels", &opt_arg))
        {
          int c = atoi (opt_arg);
          Params::raw_input_format.set_channels (c);
          Params::raw_output_format.set_channels (c);
        }
      else if (check_arg (argc, argv, &i, "--raw-rate", &opt_arg))
        {
          int r = atoi (opt_arg);
          Params::raw_input_format.set_sample_rate (r);
          Params::raw_output_format.set_sample_rate (r);
        }
      else if (check_arg (argc, argv, &i, "--quiet")
            || check_arg (argc, argv, &i, "-q"))
        {
          set_log_level (Log::WARNING);
        }
    }

  /* resort argc/argv */
  e = 1;
  for (i = 1; i < argc; i++)
    if (argv[i])
      {
        argv[e++] = argv[i];
        if (i >= e)
          argv[i] = nullptr;
      }
  *argc_p = e;
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

int
main (int argc, char **argv)
{
  parse_options (&argc, &argv);

  if (Params::have_key > 1)
    {
      error ("audiowmark: watermark key can at most be set once (--key / --test-key option)\n");
      return 1;
    }
  string op = (argc >= 2) ? argv[1] : "";

  if (op == "add" && argc == 5)
    {
      return add_watermark (argv[2], argv[3], argv[4]);
    }
  else if (op == "get" && argc == 3)
    {
      return get_watermark (argv[2], /* no ber */ "");
    }
  else if (op == "cmp" && argc == 4)
    {
      return get_watermark (argv[2], argv[3]);
    }
  else if (op == "gentest" && argc == 4)
    {
      return gentest (argv[2], argv[3]);
    }
  else if (op == "cut-start" && argc == 5)
    {
      cut_start (argv[2], argv[3], argv[4]);
    }
  else if (op == "test-subtract" && argc == 5)
    {
      test_subtract (argv[2], argv[3], argv[4]);
    }
  else if (op == "gen-key" && argc == 3)
    {
      return gen_key (argv[2]);
    }
  else
    {
      error ("audiowmark: error parsing commandline args (use audiowmark -h)\n");
      return 1;
    }
}
