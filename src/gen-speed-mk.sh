#!/bin/bash

SEEDS=$(seq 10)
STRENGTHS="10 15"
CLIPS="15 30"
MODES="0 1 2"

echo -n "all:"

for SEED in $SEEDS
do
  for MODE in $MODES
  do
    for STRENGTH in $STRENGTHS
    do
      # clips
      for CLIP in $CLIPS
      do
        echo -n " speed-$CLIP-$STRENGTH-$MODE-$SEED"
      done
      # full file
      echo -n " speed-full-$STRENGTH-$MODE-$SEED"
    done
  done
done

echo
echo

for SEED in $SEEDS
do
  for MODE in $MODES
  do
    if [ "x$MODE" == "x0" ]; then
      MODE_ARGS=""
    elif [ "x$MODE" == "x1" ]; then
      MODE_ARGS="AWM_SPEED_PATIENT=1"
    elif [ "x$MODE" == "x2" ]; then
      MODE_ARGS="AWM_TRY_SPEED=1"
    fi
    for STRENGTH in $STRENGTHS
    do
      # clips
      for CLIP in $CLIPS
      do
        FILE="speed-$CLIP-$STRENGTH-$MODE-$SEED"
        echo "$FILE:"
        echo -e "\t( cd ..; AWM_REPORT=ferv AWM_RAND_PATTERN=1 AWM_SET=huge2 AWM_PARAMS_ADD='--strength $STRENGTH' AWM_SPEED=1 $MODE_ARGS AWM_SPEED_PRE_MP3=128 AWM_CLIP='$CLIP' AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh mp3 128 ) >x$FILE"
        echo -e "\tmv x$FILE $FILE"
        echo
      done
      # full file
      FILE="speed-full-$STRENGTH-$MODE-$SEED"
      echo "$FILE:"
      echo -e "\t( cd ..; AWM_REPORT=ferv AWM_RAND_PATTERN=1 AWM_SET=huge2 AWM_PARAMS_ADD='--strength $STRENGTH' AWM_SPEED=1 $MODE_ARGS AWM_SPEED_PRE_MP3=128 AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh mp3 128 ) >x$FILE"
      echo -e "\tmv x$FILE $FILE"
      echo
    done
  done
done
