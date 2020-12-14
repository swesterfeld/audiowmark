#!/bin/bash

STRENGTHS="10 15"
CLIPS="15 30"

echo ".watermarking-speed"
echo '[frame="topbot",options="header",cols="<1,3*<"]'
echo '|=========================='

echo -n "| Strength "
for CLIP in $CLIPS
do
  echo -n "| 0:$CLIP"
done
echo -n "| 2:45"
echo

for STRENGTH in $STRENGTHS
do
  echo -n "| $STRENGTH "
  for CLIP in $CLIPS
  do
    cat speed-$CLIP-$STRENGTH-* | awk '{bad += $1; n += $2} END {if (n==0) n=1;fer=100.0*bad/n; bold=fer>0?"*":" ";printf ("| %s%.2f%s", bold, fer, bold)}'
  done
  for FULL in speed-$STRENGTH-*
  do
    if [ "$(echo $FULL | tr -d a-z0-9)" == "--" ]; then
      cat $FULL
    fi
  done | awk '{bad += $1; n += $2} END {if (n==0) n=1;fer=100.0*bad/n; bold=fer>0?"*":" ";printf ("| %s%.2f%s", bold, fer, bold)}'
  echo
done
echo '|=========================='
