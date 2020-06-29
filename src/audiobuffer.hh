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

#ifndef AUDIOWMARK_AUDIO_BUFFER_HH
#define AUDIOWMARK_AUDIO_BUFFER_HH

#include <assert.h>

class AudioBuffer
{
  const int           n_channels = 0;
  std::vector<float>  buffer;

public:
  AudioBuffer (int n_channels) :
    n_channels (n_channels)
  {
  }
  void
  write_frames (const std::vector<float>& samples)
  {
    buffer.insert (buffer.end(), samples.begin(), samples.end());
  }
  std::vector<float>
  read_frames (size_t frames)
  {
    assert (frames * n_channels <= buffer.size());
    const auto begin = buffer.begin();
    const auto end   = begin + frames * n_channels;
    std::vector<float> result (begin, end);
    buffer.erase (begin, end);
    return result;
  }
  size_t
  can_read_frames() const
  {
    return buffer.size() / n_channels;
  }
};

#endif /* AUDIOWMARK_AUDIO_BUFFER_HH */
