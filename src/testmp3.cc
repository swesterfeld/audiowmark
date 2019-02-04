#include "mp3.hh"

int
main (int argc, char **argv)
{
  WavData wd;
  mp3_try_load (argv[1], wd);
}
