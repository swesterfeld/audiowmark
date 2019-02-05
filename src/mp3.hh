#ifndef AUDIOWMARK_MP3_HH
#define AUDIOWMARK_MP3_HH

#include <string>

#include "wavdata.hh"

bool mp3_try_load (const std::string& filename, WavData& wav_data);

#endif /* AUDIOWMARK_MP3_HH */
