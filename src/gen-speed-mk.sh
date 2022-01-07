#!/bin/bash

SEEDS=$(seq 10)
STRENGTHS="10 15"
CLIPS="15 30"

echo -n "all:"

for SEED in $SEEDS
do
  for STRENGTH in $STRENGTHS
  do
    # clips
    for CLIP in $CLIPS
    do
      echo -n " speed-$CLIP-$STRENGTH-$SEED"
    done
    # full file
    echo -n " speed-full-$STRENGTH-$SEED"
  done
done

echo
echo

for SEED in $SEEDS
do
  for STRENGTH in $STRENGTHS
  do
    # clips
    for CLIP in $CLIPS
    do
      FILE="speed-$CLIP-$STRENGTH-$SEED"
      echo "$FILE:"
      echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_SET=huge2 AWM_PARAMS='--strength $STRENGTH' AWM_SPEED=1 AWM_SPEED_PRE_MP3=128 AWM_CLIP='$CLIP' AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh mp3 128 ) >x$FILE"
      echo -e "\tmv x$FILE $FILE"
      echo
    done
    # full file
    FILE="speed-full-$STRENGTH-$SEED"
    echo "$FILE:"
    echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_SET=huge2 AWM_PARAMS='--strength $STRENGTH' AWM_SPEED=1 AWM_SPEED_PRE_MP3=128 AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh mp3 128 ) >x$FILE"
    echo -e "\tmv x$FILE $FILE"
    echo
  done
done
