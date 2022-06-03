#!/bin/bash

source test-common.sh

IN_WAV=clip-decoder-test.wav
OUT_WAV=clip-decoder-test-out.wav
CUT_WAV=clip-decoder-test-out-cut.wav

audiowmark test-gen-noise $IN_WAV 30 44100
audiowmark_add $IN_WAV $OUT_WAV $TEST_MSG
audiowmark_cmp --expect-matches 1 $OUT_WAV $TEST_MSG
# cut 1 second 300 samples
audiowmark cut-start $OUT_WAV $CUT_WAV 44300
audiowmark_cmp --expect-matches 1 $CUT_WAV $TEST_MSG

rm $IN_WAV $OUT_WAV $CUT_WAV
exit 0
