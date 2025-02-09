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

class WavChunkLoader
{
  std::string                       m_filename;
  double                            m_chunk_size;
  double                            m_rate;
  double                            m_time_offset;
  std::unique_ptr<AudioInputStream> m_in_stream;
  std::vector<float>                m_samples1;
  std::vector<float>                m_samples2;
  WavData                           m_wav_data;
public:
  WavChunkLoader (const std::string& filename, double chunk_size, double rate);

  Error           load_next_chunk();
  bool            done();
  const WavData&  wav_data();
  double          time_offset();
};

#endif
