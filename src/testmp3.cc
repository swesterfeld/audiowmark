#include "mp3.hh"

int
main (int argc, char **argv)
{
  WavData wd;
  if (argc >= 2)
    {
      if (mp3_try_load (argv[1], wd))
        {
          int sec = wd.n_values() / wd.n_channels() / wd.sample_rate();

          printf ("loaded mp3 %s: %d:%02d\n", argv[1], sec / 60, sec % 60);
          if (argc == 3)
            {
              wd.save (argv[2]);
              printf ("saved wav: %s\n", argv[2]);
            }
        }
      else
        {
          printf ("mp3 try load %s failed\n", argv[1]);
          return 1;
        }
    }
}
