#!/bin/bash

SEEDS=$(seq 5)
STRENGTHS="10 15"
QUALITIES="128 256"
CLIPS="15 30 1000"

echo -n "all:"

for SEED in $SEEDS
do
  for STRENGTH in $STRENGTHS
  do
    for QUALITY in $QUALITIES
    do
      for CLIP in $CLIPS
      do
        echo -n " speed-$CLIP-$STRENGTH-$QUALITY-$SEED"
      done
    done
  done
done

echo
echo

for SEED in $SEEDS
do
  for STRENGTH in $STRENGTHS
  do
    for QUALITY in $QUALITIES
    do
      for CLIP in $CLIPS
      do
        FILE="speed-$CLIP-$STRENGTH-$QUALITY-$SEED"
        echo "$FILE:"
        echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_SET=huge2 AWM_PARAMS='--strength $STRENGTH' AWM_SPEED=1 AWM_CLIP='$CLIP' AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh mp3 $QUALITY ) >x$FILE"
        echo -e "\tmv x$FILE $FILE"
        echo
      done
    done
  done
done
