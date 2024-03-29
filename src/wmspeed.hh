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

#ifndef AUDIOWMARK_WM_SPEED_HH
#define AUDIOWMARK_WM_SPEED_HH

#include "wavdata.hh"
#include "random.hh"

struct DetectSpeedResult
{
  Key key;
  double speed = 0;
};

std::vector<DetectSpeedResult> detect_speed (const std::vector<Key>& key_list, const WavData& in_data, bool print_results);

#endif /* AUDIOWMARK_WM_SPEED_HH */
