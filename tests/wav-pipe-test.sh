#!/bin/bash

source test-common.sh

IN_WAV=wav-pipe-test.wav
OUT1_WAV=wav-pipe-test-out1.wav
OUT2_WAV=wav-pipe-test-out2.wav
OUT3_WAV=wav-pipe-test-out3.wav

for BITS in 16 24 32
do
  audiowmark test-gen-noise --bits $BITS $IN_WAV 200 44100

  [ "x$BITS" == "x$(audiowmark test-info $IN_WAV bit_depth)" ] || die "generated input bit depth is not correct"

  cat $IN_WAV   | audiowmark_add --test-key 1 --test-no-limiter --format wav-pipe - - $TEST_MSG > $OUT1_WAV || die "watermark from pipe failed"
  cat $OUT1_WAV | audiowmark_add --test-key 2 --test-no-limiter --format wav-pipe - - $TEST_MSG > $OUT2_WAV || die "watermark from pipe failed"
  cat $OUT2_WAV | audiowmark_add --test-key 3 --test-no-limiter --format wav-pipe - - $TEST_MSG > $OUT3_WAV || die "watermark from pipe failed"

  check_length $IN_WAV $OUT1_WAV
  check_length $IN_WAV $OUT2_WAV
  check_length $IN_WAV $OUT3_WAV
  check_snr $IN_WAV $OUT1_WAV 32
  check_snr $IN_WAV $OUT2_WAV 29
  check_snr $IN_WAV $OUT3_WAV 27

  audiowmark_cmp --expect-matches 0 $OUT3_WAV $TEST_MSG
  audiowmark_cmp --expect-matches 5 --test-key 1 $OUT3_WAV $TEST_MSG
  audiowmark_cmp --expect-matches 5 --test-key 2 $OUT3_WAV $TEST_MSG
  audiowmark_cmp --expect-matches 5 --test-key 3 $OUT3_WAV $TEST_MSG

  # for wav-pipe format: 16 / 24 / 32 bit input should produce the same number of output bits
  BTEST=$BITS:$(audiowmark test-info $OUT3_WAV bit_depth)
  [[ "$BTEST" =~ ^(16:16|24:24|32:32)$ ]] || die "unexpected input/output bit depth $BTEST"

  rm $IN_WAV $OUT1_WAV $OUT2_WAV $OUT3_WAV
done

exit 0
