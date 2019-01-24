#!/bin/bash

TRANSFORM=$1
if [ "x$AWM_SET" == "x" ]; then
  AWM_SET=small
fi
if [ "x$AWM_SEEDS" == "x" ]; then
  AWM_SEEDS=0
fi
if [ "x$AWM_REPORT" == "x" ]; then
  AWM_REPORT=fer
fi
if [ "x$AWM_FILE" == "x" ]; then
  AWM_FILE=t
fi

{
  if [ "x$AWM_SET" == "xsmall" ]; then
    ls test/T*
  elif [ "x$AWM_SET" == "xbig" ]; then
    cat test_list
  elif [ "x$AWM_SET" == "xhuge" ]; then
    ls huge/T*
  else
    echo "bad AWM_SET $AWM_SET" >&2
    exit 1
  fi
} | while read i
do
  for SEED in $AWM_SEEDS
  do
    echo $i

    if [ "x$AWM_RAND_PATTERN" != "x" ]; then
      # random pattern, 128 bit
      PATTERN=$(
        for i in $(seq 16)
        do
          printf "%02x" $((RANDOM % 256))
        done
      )
    else
      # pseudo random pattern, 128 bit
      PATTERN=4e1243bd22c66e76c2ba9eddc1f91394
    fi

    audiowmark add "$i" ${AWM_FILE}.wav $PATTERN $AWM_PARAMS --test-key $SEED >/dev/null
    if [ "x$AWM_RAND_CUT" != x ]; then
      CUT=$RANDOM
      audiowmark cut-start "${AWM_FILE}.wav" "${AWM_FILE}.wav" $CUT
      TEST_CUT_ARGS="--test-cut $CUT"
    else
      TEST_CUT_ARGS=""
    fi
    if [ "x$TRANSFORM" == "xmp3" ]; then
      if [ "x$2" == "x" ]; then
        echo "need mp3 bitrate" >&2
        exit 1
      fi
      lame -b $2 ${AWM_FILE}.wav ${AWM_FILE}.mp3 --quiet
      rm ${AWM_FILE}.wav
      ffmpeg -i ${AWM_FILE}.mp3 ${AWM_FILE}.wav -v quiet -nostdin

      # some (low) mpeg quality settings use a lower sample rate
      if [ "x$(soxi -r ${AWM_FILE}.wav)" != "x44100" ]; then
        sox ${AWM_FILE}.wav ${AWM_FILE}r.wav rate 44100
        mv ${AWM_FILE}r.wav ${AWM_FILE}.wav
      fi
    elif [ "x$TRANSFORM" == "xdouble-mp3" ]; then
      if [ "x$2" == "x" ]; then
        echo "need mp3 bitrate" >&2
        exit 1
      fi
      # first mp3 step (fixed bitrate)
      lame -b 128 ${AWM_FILE}.wav ${AWM_FILE}.mp3 --quiet
      rm ${AWM_FILE}.wav
      ffmpeg -i ${AWM_FILE}.mp3 ${AWM_FILE}.wav -v quiet -nostdin

      # second mp3 step
      lame -b $2 ${AWM_FILE}.wav ${AWM_FILE}.mp3 --quiet
      rm ${AWM_FILE}.wav
      ffmpeg -i ${AWM_FILE}.mp3 ${AWM_FILE}.wav -v quiet -nostdin

      # some (low) mpeg quality settings use a lower sample rate
      if [ "x$(soxi -r ${AWM_FILE}.wav)" != "x44100" ]; then
        sox ${AWM_FILE}.wav ${AWM_FILE}r.wav rate 44100
        mv ${AWM_FILE}r.wav ${AWM_FILE}.wav
      fi
    elif [ "x$TRANSFORM" == "xogg" ]; then
      if [ "x$2" == "x" ]; then
        echo "need ogg bitrate" >&2
        exit 1
      fi
      oggenc -b $2 ${AWM_FILE}.wav -o ${AWM_FILE}.ogg --quiet
      oggdec ${AWM_FILE}.ogg -o ${AWM_FILE}.wav --quiet
    elif [ "x$TRANSFORM" == "x" ]; then
      :
    else
      echo "unknown transform $TRANSFORM" >&2
      exit 1
    fi
    # blind decoding
    audiowmark cmp ${AWM_FILE}.wav $PATTERN $AWM_PARAMS --test-key $SEED $TEST_CUT_ARGS
    # decoding with original
    # audiowmark cmp-delta "$i" t.wav $PATTERN $AWM_PARAMS --test-key $SEED
  done
done | {
  if [ "x$AWM_REPORT" == "xfer" ]; then
    awk 'BEGIN { bad = n = 0 } $1 == "match_count" { if ($2 == 0) bad++; n++; } END { print bad, n, bad * 100.0 / n; }'
  elif [ "x$AWM_REPORT" == "xsync" ]; then
    awk 'BEGIN { bad = n = 0 } $1 == "sync_match" { bad += (3 - $2) / 3.0; n++; } END { print bad, n, bad * 100.0 / n; }'
  elif [ "x$AWM_REPORT" == "xsyncv" ]; then
    awk '{ print "###", $0; } $1 == "sync_match" { correct += $2; missing += 3 - $2; incorrect += $3-$2; print "correct:", correct, "missing:", missing, "incorrect:", incorrect; }'
  else
    echo "unknown report $AWM_REPORT" >&2
    exit 1
  fi
}
