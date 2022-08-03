#!/bin/bash

# set -Eeuo pipefail # -x
set -Eeo pipefail

TRANSFORM=$1
if [ "x$AWM_TRUNCATE" != "x" ]; then
  AWM_REPORT=truncv
fi
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
if [ "x$AWM_MULTI_CLIP" == "x" ]; then
  AWM_MULTI_CLIP=1
fi
if [ "x$AWM_PATTERN_BITS" == "x" ]; then
  AWM_PATTERN_BITS=128
fi

audiowmark_cmp()
{
  audiowmark cmp "$@" || {
    if [ "x$AWM_FAIL_DIR" != "x" ]; then
      mkdir -p $AWM_FAIL_DIR
      SUM=$(sha1sum $1 | awk '{print $1;}')
      cp -av $1 $AWM_FAIL_DIR/${AWM_FILE}.${SUM}.wav
    fi
  }
}

{
  if [ "x$AWM_SET" == "xsmall" ]; then
    ls test/T*
  elif [ "x$AWM_SET" == "xbig" ]; then
    cat test_list
  elif [ "x$AWM_SET" != "x" ] && [ -d "$AWM_SET" ] && [ -f "$AWM_SET/T001"*wav ]; then
    ls $AWM_SET/T*
  else
    echo "bad AWM_SET $AWM_SET" >&2
    exit 1
  fi
} | while read i
do
  for SEED in $AWM_SEEDS
  do
    echo in_file $i

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
    PATTERN=${PATTERN:0:$((AWM_PATTERN_BITS / 4))}
    echo in_pattern $PATTERN
    echo in_flags $AWM_PARAMS $AWM_PARAMS_ADD --test-key $SEED
    audiowmark add "$i" ${AWM_FILE}.wav $PATTERN $AWM_PARAMS $AWM_PARAMS_ADD --test-key $SEED --quiet
    CUT=0
    if [ "x$AWM_ALWAYS_CUT" != x ]; then
      CUT="$AWM_ALWAYS_CUT"
    fi
    if [ "x$AWM_RAND_CUT" != x ]; then
      CUT=$((CUT + RANDOM))
    fi
    if [ "x$CUT" != x0 ]; then
      audiowmark cut-start "${AWM_FILE}.wav" "${AWM_FILE}.wav" $CUT
      TEST_CUT_ARGS="--test-cut $CUT"
      echo in_cut $CUT
    else
      TEST_CUT_ARGS=""
    fi
    if [ "x$AWM_SPEED" != x ]; then
      if [ "x$AWM_SPEED_PRE_MP3" != x ]; then
        # first (optional) mp3 step: simulate quality loss before speed change
        lame -b "$AWM_SPEED_PRE_MP3" ${AWM_FILE}.wav ${AWM_FILE}.mp3 --quiet
        rm ${AWM_FILE}.wav
        ffmpeg -i ${AWM_FILE}.mp3 ${AWM_FILE}.wav -v quiet -nostdin
      fi

      [ -z $SPEED_SEED ] && SPEED_SEED=0
      SPEED=$(audiowmark test-speed $SPEED_SEED --test-key $SEED)
      SPEED_SEED=$((SPEED_SEED + 1))
      echo in_speed $SPEED

      sox -D -V1 ${AWM_FILE}.wav ${AWM_FILE}.speed.wav speed $SPEED
      mv ${AWM_FILE}.speed.wav ${AWM_FILE}.wav

      if [ "x$AWM_SPEED_PATIENT" != x ]; then
        TEST_SPEED_ARGS="--detect-speed-patient --test-speed $SPEED"
      elif [ "x$AWM_TRY_SPEED" != x ]; then
        TEST_SPEED_ARGS="--try-speed $SPEED"
      else
        TEST_SPEED_ARGS="--detect-speed --test-speed $SPEED"
      fi
    else
      TEST_SPEED_ARGS=""
    fi
    if [ "x$TRANSFORM" == "xmp3" ]; then
      if [ "x$2" == "x" ]; then
        echo "need mp3 bitrate" >&2
        exit 1
      fi
      lame -b $2 ${AWM_FILE}.wav ${AWM_FILE}.mp3 --quiet
      OUT_FILE=${AWM_FILE}.mp3
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
      OUT_FILE=${AWM_FILE}.mp3
    elif [ "x$TRANSFORM" == "xogg" ]; then
      if [ "x$2" == "x" ]; then
        echo "need ogg bitrate" >&2
        exit 1
      fi
      oggenc -b $2 ${AWM_FILE}.wav -o ${AWM_FILE}.ogg --quiet
      OUT_FILE=${AWM_FILE}.ogg
    elif [ "x$TRANSFORM" == "x" ]; then
      OUT_FILE=${AWM_FILE}.wav
    else
      echo "unknown transform $TRANSFORM" >&2
      exit 1
    fi
    echo
    if [ "x${AWM_CLIP}" != "x" ]; then
      for CLIP in $(seq $AWM_MULTI_CLIP)
      do
        audiowmark test-clip $OUT_FILE ${OUT_FILE}.clip.wav $((CLIP_SEED++)) $AWM_CLIP --test-key $SEED
        audiowmark_cmp ${OUT_FILE}.clip.wav $PATTERN $AWM_PARAMS --test-key $SEED $TEST_CUT_ARGS $TEST_SPEED_ARGS
        rm ${OUT_FILE}.clip.wav
        echo
      done
    elif [ "x$AWM_REPORT" == "xtruncv" ]; then
      for TRUNC in $AWM_TRUNCATE
      do
        audiowmark_cmp $OUT_FILE $PATTERN $AWM_PARAMS --test-key $SEED $TEST_CUT_ARGS $TEST_SPEED_ARGS --test-truncate $TRUNC | sed "s/^/$TRUNC /g"
        echo
      done
    else
      audiowmark_cmp $OUT_FILE $PATTERN $AWM_PARAMS --test-key $SEED $TEST_CUT_ARGS $TEST_SPEED_ARGS
      echo
    fi
    rm -f ${AWM_FILE}.wav $OUT_FILE # cleanup temp files
  done
done | {
  if [ "x$AWM_REPORT" == "xfer" ]; then
    awk 'BEGIN { bad = n = 0 } $1 == "match_count" { if ($2 == 0) bad++; n++; } END { print bad, n, bad * 100.0 / (n > 0 ? n : 1); }'
  elif [ "x$AWM_REPORT" == "xferv" ]; then
    awk 'BEGIN { bad = n = 0 } { print "###", $0; } $1 == "match_count" { if ($2 == 0) bad++; n++; } END { print bad, n, bad * 100.0 / (n > 0 ? n : 1); }'
  elif [ "x$AWM_REPORT" == "xsync" ]; then
    awk 'BEGIN { bad = n = 0 } $1 == "sync_match" { bad += (3 - $2) / 3.0; n++; } END { print bad, n, bad * 100.0 / (n > 0 ? n : 1); }'
  elif [ "x$AWM_REPORT" == "xsyncv" ]; then
    awk '{ print "###", $0; } $1 == "sync_match" { correct += $2; missing += 3 - $2; incorrect += $3-$2; print "correct:", correct, "missing:", missing, "incorrect:", incorrect; }'
  elif [ "x$AWM_REPORT" == "xtruncv" ]; then
    awk ' {
            print "###", $0;
          }
          $2 == "match_count" {
            if (!n[$1])
              {
                n[$1]   = 0;
                bad[$1] = 0;
              }
            if ($3 == 0)
              bad[$1]++;
            n[$1]++;
          }
          END {
            for (trunc in n) {
              print trunc, bad[trunc], n[trunc], bad[trunc] * 100.0 / n[trunc];
            }
          }'
  else
    echo "unknown report $AWM_REPORT" >&2
    exit 1
  fi
}
