# pseudo random pattern

PATTERN=4e1243bd22c66e76c2ba9eddc1f91394e57f9f83
TRANSFORM=$1

peaq_print_scores()
{
  awk '
    BEGIN {
      odg = "*ERROR*";
      di  = "*ERROR*";
    }
    /Objective Difference Grade:/ {
      odg = $NF;
    }
    /Distortion Index:/ {
      di = $NF;
    }
    END {
      print odg, di
    }'
}

for i in test/T*
do
  audiowmark add "$i" t.wav $PATTERN $AWM_PARAMS >/dev/null
  audiowmark scale "$i" t_in.wav
  echo $i $(peaq t_in.wav t.wav | peaq_print_scores)
done | awk 'BEGIN {
    max=-10;
    min=10;
    avg=0;
    count=0;
  }
  {
    if ($2 > max)
      max = $2;
    if ($2 < min)
      min = $2;
    avg += $2;
    count++;
  }
  END {
    avg /= count
    printf ("%7.3f %7.3f %7.3f\n", avg, min, max);
  }'
