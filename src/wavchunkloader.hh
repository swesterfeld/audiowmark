/*
 * Copyright (C) 2025 Stefan Westerfeld
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

#ifndef AUDIOWMARK_WAV_CHUNK_LOADER_HH
#define AUDIOWMARK_WAV_CHUNK_LOADER_HH

#include <string>

#include "utils.hh"
#include "wavdata.hh"
#include "resample.hh"

class WavChunkLoader
{
  std::string                       m_filename;
  double                            m_time_offset = 0;
  std::unique_ptr<AudioInputStream> m_in_stream;
  std::unique_ptr<ResamplerImpl>    m_resampler;
  bool                              m_resampler_in_eof = false;
  WavData                           m_wav_data;
  size_t                            m_wav_data_max_size = 0;
  size_t                            m_n_overlap_samples = 0;
  size_t                            m_n_total_samples = 0;

  enum class State
  {
    NEW,
    OPEN,
    LAST_CHUNK,
    DONE,
    ERROR
  };
  State                             m_state = State::NEW;

  Error           open();
  void            update_capacity (std::vector<float>& samples, size_t need_space, size_t max_size);
  Error           refill (std::vector<float>& samples, size_t n_values, size_t max_size, bool *eof);
public:
  WavChunkLoader (const std::string& filename);

  Error           load_next_chunk();
  bool            done();
  const WavData&  wav_data();
  double          time_offset();
  double          length();
};

#endif
