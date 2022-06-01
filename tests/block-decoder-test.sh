#!/bin/bash

source test-common.sh

IN_WAV=block-decoder-test.wav
OUT_WAV=block-decoder-test-out.wav

audiowmark test-gen-noise $IN_WAV 200 44100
audiowmark_add $IN_WAV $OUT_WAV $TEST_MSG
audiowmark_cmp $OUT_WAV $TEST_MSG

rm $IN_WAV $OUT_WAV
exit 0
