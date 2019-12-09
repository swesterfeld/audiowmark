#ifndef AUDIOWMARK_MP3_HH
#define AUDIOWMARK_MP3_HH

#include <string>

#include "wavdata.hh"

bool        mp3_detect (const std::string& filename);
void        mp3_init();

#endif /* AUDIOWMARK_MP3_HH */
