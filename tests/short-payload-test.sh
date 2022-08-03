#!/bin/bash

source test-common.sh

IN_WAV=short-playload-test.wav
OUT_WAV=short-playload-test-out.wav

audiowmark test-gen-noise $IN_WAV 200 44100
audiowmark_add --short 12 $IN_WAV $OUT_WAV abc
audiowmark_cmp --short 12 $OUT_WAV abc
audiowmark_add --short 16 $IN_WAV $OUT_WAV abcd
audiowmark_cmp --short 16 $OUT_WAV abcd
audiowmark_add --short 20 $IN_WAV $OUT_WAV abcde
audiowmark_cmp --short 20 $OUT_WAV abcde

rm $IN_WAV $OUT_WAV
exit 0
