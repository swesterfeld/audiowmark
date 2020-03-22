#!/bin/bash

SEEDS=$(seq 5)
STRENGTHS="10 15 20 30"
STR_CLIPS="5 10 15 20 25 30"
MAIN_CLIPS="5 10 15 20 30 40 50 60"
MULTI_CLIP=4

echo -n "all:"

for STRENGTH in $STRENGTHS
do
  FILE="snr-$STRENGTH"
  echo -n " $FILE"
done

for SEED in $SEEDS
do
  for STRENGTH in $STRENGTHS
  do
    for CLIP in $STR_CLIPS
    do
      FILE="str-$STRENGTH-mp3-$SEED-$CLIP"
      echo -n " $FILE"
    done
  done

  for CLIP in $MAIN_CLIPS
  do
    echo -n " main-mp3-256-$SEED-$CLIP main-mp3-128-$SEED-$CLIP main-ogg-128-$SEED-$CLIP main-double-mp3-128-$SEED-$CLIP"
  done
done

echo
echo

for SEED in $SEEDS
do
  for STRENGTH in $STRENGTHS
  do
    for CLIP in $STR_CLIPS
    do
      FILE="str-$STRENGTH-mp3-$SEED-$CLIP"
      echo "$FILE:"
      echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_SET=huge2 AWM_PARAMS_ADD='--strength $STRENGTH' AWM_CLIP='$CLIP' AWM_MULTI_CLIP='$MULTI_CLIP' AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh mp3 128 ) >x$FILE"
      echo -e "\tmv x$FILE $FILE"
      echo
    done
  done
done

for STRENGTH in $STRENGTHS
do
  FILE="snr-$STRENGTH"
  echo "$FILE:"
  echo -e "\t( cd ..; AWM_SET=huge2 AWM_PARAMS='--strength $STRENGTH' snr.sh ) >x$FILE"
  echo -e "\tmv x$FILE $FILE"
  echo
done

for SEED in $SEEDS
do
  for CLIP in $MAIN_CLIPS
  do
    FILE="main-mp3-256-$SEED-$CLIP"
    echo "$FILE:"
    echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_SET=huge2 AWM_CLIP='$CLIP' AWM_MULTI_CLIP='$MULTI_CLIP' AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh mp3 256 ) >x$FILE"
    echo -e "\tmv x$FILE $FILE"
    echo

    FILE="main-mp3-128-$SEED-$CLIP"
    echo "$FILE:"
    echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_SET=huge2 AWM_CLIP='$CLIP' AWM_MULTI_CLIP='$MULTI_CLIP' AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh mp3 128 ) >x$FILE"
    echo -e "\tmv x$FILE $FILE"
    echo

    FILE="main-ogg-128-$SEED-$CLIP"
    echo "$FILE:"
    echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_SET=huge2 AWM_CLIP='$CLIP' AWM_MULTI_CLIP='$MULTI_CLIP' AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh ogg 128 ) >x$FILE"
    echo -e "\tmv x$FILE $FILE"
    echo

    FILE="main-double-mp3-128-$SEED-$CLIP"
    echo "$FILE:"
    echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_SET=huge2 AWM_CLIP='$CLIP' AWM_MULTI_CLIP='$MULTI_CLIP' AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh double-mp3 128 ) >x$FILE"
    echo -e "\tmv x$FILE $FILE"
    echo
  done
done
