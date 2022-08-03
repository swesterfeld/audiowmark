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

#include "wmcommon.hh"
#include "fft.hh"
#include "convcode.hh"
#include "shortcode.hh"

using std::string;
using std::vector;
using std::complex;

int    Params::frames_per_bit  = 2;
double Params::water_delta     = 0.01;
bool   Params::mix             = true;
bool   Params::hard            = false; // hard decode bits? (soft decoding is better)
bool   Params::snr             = false; // compute/show snr while adding watermark
bool   Params::strict          = false;
bool   Params::detect_speed    = false;
bool   Params::detect_speed_patient = false;
double Params::try_speed       = -1;
double Params::test_speed      = -1;
int    Params::have_key        = 0;
size_t Params::payload_size    = 128;
bool   Params::payload_short   = false;
int    Params::test_cut        = 0; // for sync test
bool   Params::test_no_sync    = false; // disable sync
bool   Params::test_no_limiter = false; // disable limiter
int    Params::test_truncate   = 0;
int    Params::expect_matches  = -1;

Format Params::input_format     = Format::AUTO;
Format Params::output_format    = Format::AUTO;

RawFormat Params::raw_input_format;
RawFormat Params::raw_output_format;

int    Params::hls_bit_rate = 0;

string Params::json_output;
string Params::input_label;
string Params::output_label;

FFTAnalyzer::FFTAnalyzer (int n_channels) :
  m_n_channels (n_channels),
  m_fft_processor (Params::frame_size)
{
  m_window = gen_normalized_window (Params::frame_size);
}

/* safe to call from any thread */
vector<float>
FFTAnalyzer::gen_normalized_window (size_t n_values)
{
  vector<float> window (n_values);

  /* generate analysis window */
  double window_weight = 0;
  for (size_t i = 0; i < n_values; i++)
    {
      const double n_values_2 = n_values / 2.0;
      // const double win =  window_cos ((i - n_values_2) / n_values_2);
      const double win = window_hamming ((i - n_values_2) / n_values_2);
      //const double win = 1;
      window[i] = win;
      window_weight += win;
    }

  /* normalize window using window weight */
  for (size_t i = 0; i < n_values; i++)
    {
      window[i] *= 2.0 / window_weight;
    }
  return window;
}

vector<vector<complex<float>>>
FFTAnalyzer::run_fft (const vector<float>& samples, size_t start_index)
{
  assert (samples.size() >= (Params::frame_size + start_index) * m_n_channels);

  float *frame     = m_fft_processor.in();
  float *frame_fft = m_fft_processor.out();

  vector<vector<complex<float>>> fft_out;
  for (int ch = 0; ch < m_n_channels; ch++)
    {
      size_t pos = start_index * m_n_channels + ch;
      assert (pos + (Params::frame_size - 1) * m_n_channels < samples.size());

      /* deinterleave frame data and apply window */
      for (size_t x = 0; x < Params::frame_size; x++)
        {
          frame[x] = samples[pos] * m_window[x];
          pos += m_n_channels;
        }
      /* FFT transform */
      m_fft_processor.fft();

      /* complex<float> and frame_fft have the same layout in memory */
      const complex<float> *first = (complex<float> *) frame_fft;
      const complex<float> *last  = first + Params::frame_size / 2 + 1;
      fft_out.emplace_back (first, last);
    }

  return fft_out;
}

vector<vector<complex<float>>>
FFTAnalyzer::fft_range (const vector<float>& samples, size_t start_index, size_t frame_count)
{
  vector<vector<complex<float>>> fft_out;

  /* if there is not enough space for frame_count values, return an error (empty vector) */
  if (samples.size() < (start_index + frame_count * Params::frame_size) * m_n_channels)
    return fft_out;

  for (size_t f = 0; f < frame_count; f++)
    {
      const size_t frame_start = (f * Params::frame_size) + start_index;

      vector<vector<complex<float>>> frame_result = run_fft (samples, frame_start);
      for (auto& fr : frame_result)
        fft_out.emplace_back (std::move (fr));
    }
  return fft_out;
}

int
frame_pos (int f, bool sync)
{
  static vector<int> pos_vec;

  if (pos_vec.empty())
    {
      int frame_count = mark_data_frame_count() + mark_sync_frame_count();
      for (int i = 0; i < frame_count; i++)
        pos_vec.push_back (i);

      Random random (0, Random::Stream::frame_position);
      random.shuffle (pos_vec);
    }
  if (sync)
    {
      assert (f >= 0 && size_t (f) < mark_sync_frame_count());

      return pos_vec[f];
    }
  else
    {
      assert (f >= 0 && size_t (f) < mark_data_frame_count());

      return pos_vec[f + mark_sync_frame_count()];
    }
}

int
sync_frame_pos (int f)
{
  return frame_pos (f, true);
}

int
data_frame_pos (int f)
{
  return frame_pos (f, false);
}

size_t
mark_data_frame_count()
{
  return code_size (ConvBlockType::a, Params::payload_size) * Params::frames_per_bit;
}

size_t
mark_sync_frame_count()
{
  return Params::sync_bits * Params::sync_frames_per_bit;
}

vector<MixEntry>
gen_mix_entries()
{
  const int frame_count = mark_data_frame_count();
  vector<MixEntry> mix_entries (frame_count * Params::bands_per_frame);

  UpDownGen up_down_gen (Random::Stream::data_up_down);
  int entry = 0;
  for (int f = 0; f < frame_count; f++)
    {
      const int index = data_frame_pos (f);
      UpDownArray up, down;
      up_down_gen.get (f, up, down);

      assert (up.size() == down.size());
      for (size_t i = 0; i < up.size(); i++)
        mix_entries[entry++] = { index, up[i], down[i] };
    }
  Random random (/* seed */ 0, Random::Stream::mix);
  random.shuffle (mix_entries);

  return mix_entries;
}

int
frame_count (const WavData& wav_data)
{
  return wav_data.n_values() / wav_data.n_channels() / Params::frame_size;
}

vector<int>
parse_payload (const string& bits)
{
  auto bitvec = bit_str_to_vec (bits);
  if (bitvec.empty())
    {
      error ("audiowmark: cannot parse bits '%s'\n", bits.c_str());
      return {};
    }
  if ((Params::payload_short || Params::strict) && bitvec.size() != Params::payload_size)
    {
      error ("audiowmark: number of message bits must match payload size (%zd bits)\n", Params::payload_size);
      return {};
    }
  if (bitvec.size() > Params::payload_size)
    {
      error ("audiowmark: number of bits in message '%s' larger than payload size\n", bits.c_str());
      return {};
    }
  if (bitvec.size() < Params::payload_size)
    {
      /* expand message automatically; good for testing, not so good in production (disabled by --strict) */
      vector<int> expanded_bitvec;
      for (size_t i = 0; i < Params::payload_size; i++)
        expanded_bitvec.push_back (bitvec[i % bitvec.size()]);
      bitvec = expanded_bitvec;
    }
  return bitvec;
}
