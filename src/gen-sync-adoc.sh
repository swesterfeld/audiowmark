#!/bin/bash

echo ".sync-codec-resistence"
echo '[frame="topbot",options="header",cols="<2,6*>1"]'
echo '|=========================='
echo -n "| "
for D in $(seq 10 -1 5)
do
  DELTA=$(printf "0.0%02d\n" $D)
  echo -n "| $DELTA"
done
echo
for TEST in mp3 double-mp3 ogg
do
  if [ $TEST == mp3 ]; then
    echo -n "| mp3 128kbit/s"
  elif [ $TEST == double-mp3 ]; then
    echo -n "| double mp3 128kbit/s"
  elif [ $TEST == ogg ]; then
    echo -n "| ogg 128kbit/s"
  else
    echo "error: bad TEST $TEST ???"
    exit 1
  fi
  for D in $(seq 10 -1 5)
  do
    DELTA=$(printf "0.0%02d\n" $D)
    cat $DELTA-$TEST-* | awk '{bad += $1; n += $2} END {if (n==0) n=1;fer=100.0*bad/n; bold=fer>0?"*":" ";printf ("| %s%.2f%s", bold, fer, bold)}'
  done
  echo
done
echo
echo '|=========================='
