#!/bin/bash

source test-common.sh

IN_WAV=detect-speed-test.wav
OUT_WAV=detect-speed-test-out.wav
OUTS_WAV=detect-speed-test-out-spd.wav

$AUDIOWMARK test-gen-noise detect-speed-test.wav 30 44100 || die "failed to generate noise"
for SPEED in 0.9764 1.01
do
  audiowmark_add $IN_WAV $OUT_WAV $TEST_MSG
  $AUDIOWMARK test-change-speed $OUT_WAV $OUTS_WAV $SPEED || die "failed to change speed"
  audiowmark_cmp $OUTS_WAV $TEST_MSG --detect-speed --test-speed $SPEED
  audiowmark_cmp $OUTS_WAV $TEST_MSG --detect-speed-patient --test-speed $SPEED
done

rm $IN_WAV $OUT_WAV $OUTS_WAV
exit 0
