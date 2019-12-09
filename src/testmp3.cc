#include "mp3.hh"
#include "mp3inputstream.hh"

using std::string;

int
main (int argc, char **argv)
{
  WavData wd;
  if (argc >= 2)
    {
      if (mp3_detect (argv[1]))
        {
          MP3InputStream m3i;
          Error err = m3i.open (argv[1]);
          if (err)
            {
              printf ("mp3 open %s failed: %s\n", argv[1], err.message());
              return 1;
            }

          err = wd.load (&m3i);
          if (!err)
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
              printf ("mp3 load %s failed: %s\n", argv[1], err.message());
              return 1;
            }
        }
      else
        {
          printf ("mp3 detect %s failed\n", argv[1]);
          return 1;
        }
    }
}
