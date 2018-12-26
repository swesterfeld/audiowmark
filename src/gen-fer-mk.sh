#!/bin/bash

DELTA_RANGE="0.005 0.006 0.007 0.008 0.009 0.010 0.011 0.012 0.013 0.014 0.015"
SEEDS="$(seq 0 19)"

echo -n "all: "
for SEED in $SEEDS
do
  for DELTA in $DELTA_RANGE
  do
    echo -n "$DELTA-ogg-$SEED $DELTA-mp3-$SEED $DELTA-double-mp3-$SEED "
  done
done

echo
echo

for SEED in $SEEDS
do
  for DELTA in $DELTA_RANGE
  do
    FILE="$DELTA-ogg-$SEED"
    echo "$FILE:"
    echo -e "\t( cd ..; AWM_SET=huge AWM_PARAMS='--water-delta $DELTA' AWM_REPORT=fer AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-ogg.sh ) >x$FILE"
    echo -e "\tmv x$FILE $FILE"
    echo

    FILE="$DELTA-mp3-$SEED"
    echo "$FILE:"
    echo -e "\t( cd ..; AWM_SET=huge AWM_PARAMS='--water-delta $DELTA' AWM_REPORT=fer AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-mp3.sh ) >x$FILE"
    echo -e "\tmv x$FILE $FILE"
    echo

    FILE="$DELTA-double-mp3-$SEED"
    echo "$FILE:"
    echo -e "\t( cd ..; AWM_SET=huge AWM_PARAMS='--water-delta $DELTA' AWM_REPORT=fer AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-double-mp3.sh ) >x$FILE"
    echo -e "\tmv x$FILE $FILE"
    echo
  done
done
