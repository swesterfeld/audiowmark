#!/bin/bash

STRENGTHS="10 15 20 30"
STR_CLIPS="5 10 15 20 25 30"
MAIN_CLIPS="5 10 15 20 30 40 50 60"

echo ".sync-codec-resistence$TRUNC"
echo '[frame="topbot",options="header",cols="<2,6*>1"]'
echo '|=========================='
echo -n "| "
for STRENGTH in $STRENGTHS
do
  echo -n "| $STRENGTH"
done
echo
for CLIP in $STR_CLIPS
do
  for STRENGTH in $STR_CLIPS
  do
    cat str-$STRENGTH-mp3-*-$CLIP | awk '{bad += $2; n += $3} END {if (n==0) n=1;fer=100.0*bad/n; bold=fer>0?"*":" ";printf ("| %s%.2f%s", bold, fer, bold)}'
  done
  echo
done
echo
echo '|=========================='
