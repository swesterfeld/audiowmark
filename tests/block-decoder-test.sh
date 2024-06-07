#!/bin/bash

source test-common.sh

IN_WAV=block-decoder-test.wav
OUT_WAV=block-decoder-test-out.wav

audiowmark test-gen-noise $IN_WAV 200 44100
audiowmark_add $IN_WAV $OUT_WAV $TEST_MSG
audiowmark_cmp --expect-matches 5 $OUT_WAV $TEST_MSG

check_length $IN_WAV $OUT_WAV

audiowmark_add --test-no-limiter $IN_WAV $OUT_WAV $TEST_MSG
audiowmark_cmp --expect-matches 5 $OUT_WAV $TEST_MSG

check_length $IN_WAV $OUT_WAV
check_snr $IN_WAV $OUT_WAV 32.4

rm $IN_WAV $OUT_WAV
exit 0
