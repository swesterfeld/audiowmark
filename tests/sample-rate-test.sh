#!/bin/bash

source test-common.sh

IN_WAV=sample-rate-test.wav
OUT_WAV=sample-rate-test-out.wav
OUT_48000_WAV=sample-rate-test-out-48000.wav

audiowmark test-gen-noise $IN_WAV 200 32000
audiowmark_add $IN_WAV $OUT_WAV $TEST_MSG
audiowmark_cmp --expect-matches 5 $OUT_WAV $TEST_MSG
audiowmark test-resample $OUT_WAV $OUT_48000_WAV 48000
audiowmark_cmp --expect-matches 5 $OUT_48000_WAV $TEST_MSG

rm $IN_WAV $OUT_WAV $OUT_48000_WAV
exit 0
