#include "utils.hh"
#include "random.hh"

#include <sys/time.h>

using std::vector;
using std::string;

static double
gettime()
{
  timeval tv;
  gettimeofday (&tv, 0);

  return tv.tv_sec + tv.tv_usec / 1000000.0;
}

int
main (int argc, char **argv)
{
  Random rng (0xf00f1234b00b5678U, Random::Stream::bit_order);
  for (size_t i = 0; i < 20; i++)
    {
      uint64_t x = rng();
      printf ("%016lx\n", x);
    }

  uint64_t s = 0;
  double t_start = gettime();
  size_t runs = 25000000;
  for (size_t i = 0; i < runs; i++)
    {
      s += rng();
    }
  double t_end = gettime();
  printf ("s=%016lx\n\n", s);

  printf ("%f Mvalues/sec\n", runs / (t_end - t_start) / 1000000);
}
