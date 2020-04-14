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

#ifndef AUDIOWMARK_AUDIO_STREAM_HH
#define AUDIOWMARK_AUDIO_STREAM_HH

#include <vector>
#include <memory>
#include "utils.hh"

class AudioStream
{
public:
  virtual int     bit_depth()   const = 0;
  virtual int     sample_rate() const = 0;
  virtual int     n_channels()  const = 0;

  virtual ~AudioStream();
};

class AudioInputStream : public AudioStream
{
public:
  static std::unique_ptr<AudioInputStream> create (const std::string& filename, Error& err);

  // for streams that do not know the number of frames in advance (i.e. raw input stream)
  static constexpr size_t N_FRAMES_UNKNOWN = ~size_t (0);
  virtual size_t n_frames() const = 0;

  virtual Error read_frames (std::vector<float>& samples, size_t count) = 0;
};

class AudioOutputStream : public AudioStream
{
public:
  static std::unique_ptr<AudioOutputStream> create (const std::string& filename,
    int n_channels, int sample_rate, int bit_depth, size_t n_frames, Error& err);

  virtual Error write_frames (const std::vector<float>& frames) = 0;
  virtual Error close() = 0;
};

#endif /* AUDIOWMARK_AUDIO_STREAM_HH */
