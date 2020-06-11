#!/bin/bash

STRENGTHS="10 15"
STR_CLIPS="6 10 15 20 25 30"
QUALITIES="128 256"

for LS in long short
do
  echo ".watermarking-with-$LS-payload"
  echo '[frame="topbot",options="header",cols="<1,7*>1"]'
  echo '|=========================='

  echo -n "| Strength | Quality "
  for CLIP in $STR_CLIPS
  do
    echo -n "| $CLIP"
  done
  echo

  for STRENGTH in $STRENGTHS
  do
    for QUALITY in $QUALITIES
    do
      echo -n "| $STRENGTH | $QUALITY "
      for CLIP in $STR_CLIPS
      do
        cat $LS-$CLIP-$STRENGTH-$QUALITY-* | awk '{bad += $1; n += $2} END {if (n==0) n=1;fer=100.0*bad/n; bold=fer>0?"*":" ";printf ("| %s%.2f%s", bold, fer, bold)}'
      done
      echo
    done
  done
  echo '|=========================='
done
