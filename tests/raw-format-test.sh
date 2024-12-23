source test-common.sh

IN_WAV=test-raw-format.wav
OUT_WAV=test-raw-format-out.wav
OUT2_WAV=test-raw-format-out2.wav

raw_test()
{
  FFMPEG_FMT="$1"
  shift
  AWM_FMT="$@"

  rm -f $IN_WAV $OUT_WAV $OUT2_WAV

  audiowmark test-gen-noise --bits 32 $IN_WAV 200 44100

  ffmpeg -v quiet -i $IN_WAV -f $FFMPEG_FMT -c:a pcm_$FFMPEG_FMT - | \
    audiowmark_add - - $TEST_MSG --format raw --raw-rate 44100 $AWM_FMT --test-no-limiter | \
    ffmpeg -v quiet -f $FFMPEG_FMT -ar 44100 -ac 2 -i - $OUT_WAV

  audiowmark_cmp --expect-matches 5 $OUT_WAV $TEST_MSG
  check_snr $IN_WAV $OUT_WAV 32.4

  ffmpeg -v quiet -i $IN_WAV -f $FFMPEG_FMT -c:a pcm_$FFMPEG_FMT - | \
    audiowmark_add - $OUT2_WAV $TEST_MSG --input-format raw --raw-rate 44100 $AWM_FMT --test-no-limiter

  check_length $IN_WAV $OUT_WAV
  check_length $IN_WAV $OUT2_WAV

  rm -f $IN_WAV $OUT_WAV $OUT2_WAV
}

raw_test s16le
raw_test s24le --raw-bits 24
raw_test s32le --raw-bits 32
#raw_test u16le --raw-encoding unsigned
#raw_test u24le --raw-bits 24 --raw-encoding unsigned
#raw_test u32le --raw-bits 32 --raw-encoding unsigned
raw_test f32le --raw-encoding float
raw_test f64le --raw-encoding double
raw_test s16be --raw-endian big
raw_test s24be --raw-bits 24 --raw-endian big
raw_test s32be --raw-bits 32 --raw-endian big
#raw_test u16be --raw-encoding unsigned --endian big
#raw_test u24be --raw-bits 24 --raw-encoding unsigned --endian big
#raw_test u32be --raw-bits 32 --raw-encoding unsigned --endian big
raw_test f32be --raw-encoding float --raw-endian big
raw_test f64be --raw-encoding double --raw-endian big
