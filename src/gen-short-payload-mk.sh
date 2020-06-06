#!/bin/bash

SEEDS=$(seq 5)
STRENGTHS="10 15"
QUALITIES="128 256"
CLIPS="6 10 15 20 25 30"
MULTI_CLIP=4

echo -n "all:"

for SEED in $SEEDS
do
  for STRENGTH in $STRENGTHS
  do
    for QUALITY in $QUALITIES
    do
      for CLIP in $CLIPS
      do
        echo -n " long-$CLIP-$STRENGTH-$QUALITY-$SEED short-$CLIP-$STRENGTH-$QUALITY-$SEED"
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
        FILE="long-$CLIP-$STRENGTH-$QUALITY-$SEED"
        echo "$FILE:"
        echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_ALWAYS_CUT=500000 AWM_SET=huge2 AWM_PARAMS='--strength $STRENGTH' AWM_CLIP='$CLIP' AWM_MULTI_CLIP='$MULTI_CLIP' AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh mp3 $QUALITY ) >x$FILE"
        echo -e "\tmv x$FILE $FILE"
        echo

        FILE="short-$CLIP-$STRENGTH-$QUALITY-$SEED"
        echo "$FILE:"
        echo -e "\t( cd ..; AWM_PATTERN_BITS=12 AWM_RAND_PATTERN=1 AWM_ALWAYS_CUT=500000 AWM_SET=huge2 AWM_PARAMS='--short --strength $STRENGTH' AWM_CLIP='$CLIP' AWM_MULTI_CLIP='$MULTI_CLIP' AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh mp3 $QUALITY ) >x$FILE"
        echo -e "\tmv x$FILE $FILE"
        echo
      done
    done
  done
done
