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

#ifndef AUDIOWMARK_HLS_HH
#define AUDIOWMARK_HLS_HH

#include <string>

int hls_add (const std::string& infile, const std::string& outfile, const std::string& bits);
int hls_prepare (const std::string& in_dir, const std::string& out_dir, const std::string& filename, const std::string& audio_master);

Error ff_decode (const std::string& filename, WavData& out_wav_data);

#endif /* AUDIOWMARK_MPEGTS_HH */
