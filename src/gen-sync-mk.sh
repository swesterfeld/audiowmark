#!/bin/bash

DELTA_RANGE="0.005 0.006 0.007 0.008 0.009 0.010"
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
    echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_RAND_CUT=1 AWM_SET=huge2 AWM_PARAMS='--water-delta $DELTA' AWM_REPORT=fer AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh ogg 128 ) >x$FILE"
    echo -e "\tmv x$FILE $FILE"
    echo

    FILE="$DELTA-mp3-$SEED"
    echo "$FILE:"
    echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_RAND_CUT=1 AWM_SET=huge2 AWM_PARAMS='--water-delta $DELTA' AWM_REPORT=fer AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh mp3 128 ) >x$FILE"
    echo -e "\tmv x$FILE $FILE"
    echo

    FILE="$DELTA-double-mp3-$SEED"
    echo "$FILE:"
    echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_RAND_CUT=1 AWM_SET=huge2 AWM_PARAMS='--water-delta $DELTA' AWM_REPORT=fer AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh double-mp3 128 ) >x$FILE"
    echo -e "\tmv x$FILE $FILE"
    echo
  done
done
