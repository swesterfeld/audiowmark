#!/bin/bash

source test-common.sh

IN_WAV=pipe-test.wav
OUT_WAV=pipe-test-out.wav

audiowmark test-gen-noise $IN_WAV 200 44100
cat $IN_WAV | audiowmark_add - - $TEST_MSG > $OUT_WAV || die "watermark from pipe failed"
audiowmark_cmp $OUT_WAV $TEST_MSG
cat $OUT_WAV | audiowmark_cmp - $TEST_MSG || die "watermark detection from pipe failed"

rm $IN_WAV $OUT_WAV
exit 0
