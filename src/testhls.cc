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
#include "hls.hh"
#include "sfinputstream.hh"
#include "hlsoutputstream.hh"

using std::string;
using std::regex;
using std::vector;
using std::map;
using std::min;

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
  write_frames (const std::vector<float>& frames) override
  {
    samples.insert (samples.end(), frames.begin(), frames.end());
    return Error::Code::NONE;
  }
  Error
  close() override
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
seek_perf (int sample_rate, double seconds)
{
  vector<float> samples (100);
  WavData wav_data (samples, 2, sample_rate, 16);

  double start_time = get_time();

  int rc = mark_zexpand (wav_data, seconds * sample_rate, "0c");
  if (rc != 0)
    return rc;

  double end_time = get_time();

  info ("\n\n");
  info ("total time %7.3f sec\n", end_time - start_time);
  info ("per second %7.3f ms\n", (end_time - start_time) / seconds * 1000);

  return 0;
}

int
main (int argc, char **argv)
{
  if (argc == 6 && strcmp (argv[1], "test-seek") == 0)
    {
      return test_seek (argv[2], argv[3], atoi (argv[4]), argv[5]);
    }
  else if (argc == 4 && strcmp (argv[1], "seek-perf") == 0)
    {
      return seek_perf (atoi (argv[2]), atof (argv[3]));
    }
  else if (argc == 4 && strcmp (argv[1], "ff-decode") == 0)
    {
      WavData wd;
      Error err = ff_decode (argv[2], wd);
      if (err)
        {
          error ("audiowmark: hls: ff_decode failed: %s\n", err.message());
          return 1;
        }
      err = wd.save (argv[3]);
      if (err)
        {
          error ("audiowmark: hls: save failed: %s\n", err.message());
          return 1;
        }
      return 0;
    }
  else
    {
      error ("testhls: error parsing command line arguments\n");
      return 1;
    }
}

