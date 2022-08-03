#!/bin/bash

source test-common.sh

IN_WAV=key-test.wav
KEY1=key-test-1.key
KEY2=key-test-2.key
OUT1_WAV=key-test-out1.wav
OUT2_WAV=key-test-out2.wav

TEST_MSG2=0123456789abcdef0123456789abcdef

audiowmark test-gen-noise $IN_WAV 30 44100
audiowmark gen-key $KEY1
audiowmark gen-key $KEY2

audiowmark_add --key $KEY1 $IN_WAV $OUT1_WAV $TEST_MSG
audiowmark_add --key $KEY2 $IN_WAV $OUT2_WAV $TEST_MSG2

# shouldn't be able to detect watermark without correct key

audiowmark_cmp --key $KEY1 --expect-matches 1 $OUT1_WAV $TEST_MSG
audiowmark_cmp --key $KEY2 --expect-matches 0 $OUT1_WAV $TEST_MSG
audiowmark_cmp --expect-matches 0 $OUT1_WAV $TEST_MSG

audiowmark_cmp --key $KEY2 --expect-matches 1 $OUT2_WAV $TEST_MSG2
audiowmark_cmp --key $KEY1 --expect-matches 0 $OUT2_WAV $TEST_MSG2
audiowmark_cmp --expect-matches 0 $OUT2_WAV $TEST_MSG2

rm $OUT1_WAV $OUT2_WAV

# double watermark with two different keys
audiowmark_add $IN_WAV $OUT1_WAV $TEST_MSG
audiowmark_add --test-key 42 $OUT1_WAV $OUT2_WAV $TEST_MSG2

audiowmark_cmp --expect-matches 1 $OUT2_WAV $TEST_MSG
audiowmark_cmp --test-key 42 --expect-matches 1 $OUT2_WAV $TEST_MSG2

rm $IN_WAV $KEY1 $KEY2 $OUT1_WAV $OUT2_WAV
exit 0
