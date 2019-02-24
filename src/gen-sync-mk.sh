#!/bin/bash

STRENGTH_RANGE=$(seq 5 10)
SEEDS="$(seq 0 19)"

echo -n "all: "
for SEED in $SEEDS
do
  for STRENGTH in $STRENGTH_RANGE
  do
    echo -n "$STRENGTH-ogg-$SEED $STRENGTH-mp3-$SEED $STRENGTH-double-mp3-$SEED "
  done
done

echo
echo

for SEED in $SEEDS
do
  for STRENGTH in $STRENGTH_RANGE
  do
    FILE="$STRENGTH-ogg-$SEED"
    echo "$FILE:"
    echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_RAND_CUT=1 AWM_SET=huge2 AWM_PARAMS='--strength $STRENGTH' AWM_REPORT=ferv AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh ogg 128 ) >x$FILE"
    echo -e "\tmv x$FILE $FILE"
    echo

    FILE="$STRENGTH-mp3-$SEED"
    echo "$FILE:"
    echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_RAND_CUT=1 AWM_SET=huge2 AWM_PARAMS='--strength $STRENGTH' AWM_REPORT=ferv AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh mp3 128 ) >x$FILE"
    echo -e "\tmv x$FILE $FILE"
    echo

    FILE="$STRENGTH-double-mp3-$SEED"
    echo "$FILE:"
    echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_RAND_CUT=1 AWM_SET=huge2 AWM_PARAMS='--strength $STRENGTH' AWM_REPORT=ferv AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh double-mp3 128 ) >x$FILE"
    echo -e "\tmv x$FILE $FILE"
    echo
  done
done
