#!/bin/bash

source test-common.sh

IN_WAV=clip-decoder-test.wav
OUT_WAV=clip-decoder-test-out.wav
CUT_WAV=clip-decoder-test-out-cut.wav

$AUDIOWMARK test-gen-noise $IN_WAV 30 44100 || die "failed to generate noise"
audiowmark_add $IN_WAV $OUT_WAV $TEST_MSG
audiowmark_cmp $OUT_WAV $TEST_MSG
$AUDIOWMARK cut-start $OUT_WAV $CUT_WAV 44300 || die "failed to cut input"
audiowmark_cmp $CUT_WAV $TEST_MSG

rm $IN_WAV $OUT_WAV $CUT_WAV
exit 0
