#!/bin/bash


source test-common.sh

TESTWAVFORMAT=$TOP_BUILDDIR/src/testwavformat

testwavformat()
{
  if [ "x$Q" == "x1" ] && [ -z "$V" ]; then
    : # silent
  else
    echo >&2 ==== testwavformat "$@" ====
  fi
  $TESTWAVFORMAT "$@" || die "failed to run testwavformat $@"
}

compare_fmt_snr()
{
  INFMT=$(testwavformat detect $1)
  OUTFMT=$(testwavformat detect $2)
  EXPECTFMT=$(echo $INFMT | sed s/pcm_8/pcm_16/g)
  if [ "x$Q" == "x1" ] && [ -z "$V" ]; then
    : # silent
  else
    echo >&2 ==== infmt $INFMT outfmt $OUTFMT expectfmt $EXPECTFMT ====
  fi
  if [ "x$EXPECTFMT" != "x$OUTFMT" ]; then
    die "format mismatch $EXPECTFMT $OUTFMT"
    exit 1
  fi
  check_snr $1 $2 $3
}

IN_WAV=wav-subformat-in.wav
FMT_WAV=wav-subformat.wav
MARK_WAV=wav-subformat-mark.wav

audiowmark test-gen-noise --bits 32 $IN_WAV 200 44100

for FMT in $(testwavformat list)
do
  testwavformat convert $IN_WAV $FMT_WAV $FMT

  # FILE -> FILE
  audiowmark_add $FMT_WAV $MARK_WAV $TEST_MSG --test-no-limiter
  compare_fmt_snr $FMT_WAV $MARK_WAV 32.3

  # STDIN -> FILE
  cat $FMT_WAV | audiowmark_add - $MARK_WAV $TEST_MSG --test-no-limiter
  compare_fmt_snr $FMT_WAV $MARK_WAV 32.3

  # FILE -> STDOUT
  audiowmark_add $FMT_WAV - $TEST_MSG --test-no-limiter > $MARK_WAV
  compare_fmt_snr $FMT_WAV $MARK_WAV 32.3

  # STDIN (WavPipe) -> FILE
  cat $FMT_WAV | audiowmark_add - $MARK_WAV $TEST_MSG --test-no-limiter --format wav-pipe
  compare_fmt_snr $FMT_WAV $MARK_WAV 32.3

  # STDIN (WavPipe) -> STDOUT (WavPipe)
  cat $FMT_WAV | audiowmark_add - - $TEST_MSG --test-no-limiter --format wav-pipe > $MARK_WAV
  compare_fmt_snr $FMT_WAV $MARK_WAV 32.3
done

rm $IN_WAV $FMT_WAV $MARK_WAV
