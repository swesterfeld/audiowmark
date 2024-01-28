#!/bin/bash

SEEDS=$(seq 5)
STRENGTHS="10 15"
MP3_QUALITIES="128 64 48"
MULTI_CLIP=4

if [ "x$1" == "xmk" ]; then
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
      for Q in $MP3_QUALITIES
      do
        FILE="mp3-$Q-$STRENGTH-$SEED"
        echo -n " $FILE"
      done
    done
  done

  echo
  echo

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
    for STRENGTH in $STRENGTHS
    do
      for Q in $MP3_QUALITIES
      do
        FILE="mp3-$Q-$STRENGTH-$SEED"
        echo "$FILE:"
        echo -e "\t( cd ..; AWM_RAND_PATTERN=1 AWM_SET=huge2 AWM_CLIP='50' AWM_PARAMS_ADD='--strength $STRENGTH' AWM_MULTI_CLIP='$MULTI_CLIP' AWM_SEEDS=$SEED AWM_FILE='t-$FILE' ber-test.sh mp3 $Q ) >x$FILE"
        echo -e "\tmv x$FILE $FILE"
        echo
      done
    done
  done
fi

if [ "x$1" == "xadoc" ]; then
  echo '[frame="topbot",options="header",cols="<1,4*<"]'
  echo '|=========================='
  echo '| *Strength* | *SNR* | *mp3 128kbit/s* | *mp3 64kbit/s* | *mp3 48kbit/s*'

  for STRENGTH in $STRENGTHS
  do
    echo -n "| *$STRENGTH* "
    cat snr-$STRENGTH | awk '{printf ("| %.2f ", $1);}'
    for Q in $MP3_QUALITIES
    do
      for FILE in "mp3-$Q-$STRENGTH-*"
      do
        cat $FILE
      done | grep -v '#' | awk '{bad += $1; n += $2} END {if (n==0) n=1;fer=100.0*bad/n; bold=fer>0?"*":" ";printf ("| %s%.2f%s", bold, fer, bold)}'
    done
    echo
  done

  echo '|=========================='
fi
