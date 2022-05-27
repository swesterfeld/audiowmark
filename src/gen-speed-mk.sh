#!/bin/bash

SEEDS=$(seq 10)
STRENGTHS="10 15"
CLIPS="15 30"
PMODES="0 1"

echo -n "all:"

for SEED in $SEEDS
do
  for PMODE in $PMODES
  do
    for STRENGTH in $STRENGTHS
    do
      # clips
      for CLIP in $CLIPS
      do
        echo -n " speed-$CLIP-$STRENGTH-$PMODE-$SEED"
      done
      # full file
      echo -n " speed-full-$STRENGTH-$PMODE-$SEED"
    done
  done
done

echo
echo

for SEED in $SEEDS
do
  for PMODE in $PMODES
  do
    if [ "x$PMODE" == "x0" ]; then
      PATIENT=""
    else
      PATIENT="AWM_SPEED_PATIENT=1"
    fi
    for STRENGTH in $STRENGTHS
    do
      # clips
      for CLIP in $CLIPS
      do
        FILE="speed-$CLIP-$STRENGTH-$PMODE-$SEED"
        echo "$FILE:"
        echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_SET=huge2 AWM_PARAMS='--strength $STRENGTH' AWM_SPEED=1 $PATIENT AWM_SPEED_PRE_MP3=128 AWM_CLIP='$CLIP' AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh mp3 128 ) >x$FILE"
        echo -e "\tmv x$FILE $FILE"
        echo
      done
      # full file
      FILE="speed-full-$STRENGTH-$PMODE-$SEED"
      echo "$FILE:"
      echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_SET=huge2 AWM_PARAMS='--strength $STRENGTH' AWM_SPEED=1 $PATIENT AWM_SPEED_PRE_MP3=128 AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh mp3 128 ) >x$FILE"
      echo -e "\tmv x$FILE $FILE"
      echo
    done
  done
done
