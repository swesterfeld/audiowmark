#!/bin/bash

for TEST in mp3 double-mp3 ogg
do
  echo ".$TEST"
  echo '[frame="topbot",options="header",cols="12*>"]'
  echo '|=========================='
  echo -n "| "
  for D in $(seq 15 -1 5)
  do
    DELTA=$(printf "0.0%02d\n" $D)
    echo -n "| $DELTA"
  done
  echo
  for BITRATE in 512 256 196 128 96 64
  do
    echo -n "| $BITRATE" | sed 's/512/wav/g'
    for D in $(seq 15 -1 5)
    do
      DELTA=$(printf "0.0%02d\n" $D)
      cat $DELTA-$TEST-* | grep ^$BITRATE | awk '{bad += $2; n += $3} END {fer=100.0*bad/n; bold=fer>0?"*":" ";printf ("| %s%.2f%s", bold, fer, bold)}'
    done
    echo
  done
  echo '|=========================='
  echo
done
