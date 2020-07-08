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

#include "utils.hh"
#include "mpegts.hh"
#include "sfinputstream.hh"
#include "hlsoutputstream.hh"
#include "wmcommon.hh"

using std::string;
using std::map;
using std::min;

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
  size_t prev_ctx = min<size_t> (1024 * 3, prev_size);

  info ("hls_time_elapsed_decode %f\n", (get_time() - start_time1) * 1000 /* ms */);
  start_time1 = get_time();

  HLSOutputStream out_stream (in_stream.n_channels(), in_stream.sample_rate(), in_stream.bit_depth());

  int bit_rate = Params::hls_bit_rate ? Params::hls_bit_rate : 256000;
  out_stream.set_bit_rate (Params::hls_bit_rate ? Params::hls_bit_rate : 256000);
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


