#!/bin/bash

source test-common.sh

if [ "x$Q" == "x1" ] && [ -z "$V" ]; then
  FFMPEG_Q="-v quiet"
fi

set -e

HLS_DIR=hls-test-dir.$$
mkdir -p $HLS_DIR

# generate input sample
audiowmark test-gen-noise $HLS_DIR/test-input.wav 200 44100

# convert to hls
ffmpeg $FFMPEG_Q -i $HLS_DIR/test-input.wav \
  -f hls \
  -c:a:0 aac -ab 192k \
  -master_pl_name replay.m3u8 \
  -hls_list_size 0 -hls_time 10 $HLS_DIR/as%v/out.m3u8

# prepare hls segments for watermarking
audiowmark hls-prepare $HLS_DIR/as0 $HLS_DIR/as0prep out.m3u8 $HLS_DIR/test-input.wav

# watermark hls segments individually
mkdir -p $HLS_DIR/as0m
for i in $(cd $HLS_DIR/as0; ls out*.ts)
do
  audiowmark hls-add $HLS_DIR/as0prep/$i $HLS_DIR/as0m/$i $TEST_MSG
done
cp $HLS_DIR/as0/out.m3u8 $HLS_DIR/as0m/out.m3u8

# convert watermarked hls back to wav
ffmpeg $FFMPEG_Q -y -i $HLS_DIR/as0m/out.m3u8 $HLS_DIR/test-output.wav

# detect watermark from wav
audiowmark_cmp --expect-matches 5 $HLS_DIR/test-output.wav $TEST_MSG

rm $HLS_DIR/as0*/*.ts
rm $HLS_DIR/as0*/out.m3u8
rmdir $HLS_DIR/as0*
rm $HLS_DIR/test-*.wav
rm $HLS_DIR/replay.m3u8
rmdir $HLS_DIR

exit 0
