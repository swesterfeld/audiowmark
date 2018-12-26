# pseudo random pattern, 128 bit

PATTERN=4e1243bd22c66e76c2ba9eddc1f91394
TRANSFORM=$1
if [ "x$AWM_SET" == "x" ]; then
  AWM_SET=small
fi
if [ "x$AWM_SEEDS" == "x" ]; then
  AWM_SEEDS=0
fi
if [ "x$AWM_REPORT" == "x" ]; then
  AWM_REPORT=ber
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
    audiowmark add "$i" ${AWM_FILE}.wav $PATTERN $AWM_PARAMS --seed $SEED >/dev/null
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
    audiowmark cmp ${AWM_FILE}.wav $PATTERN $AWM_PARAMS --seed $SEED
    # decoding with original
    # audiowmark cmp-delta "$i" t.wav $PATTERN $AWM_PARAMS --seed $SEED
  done
done | grep bit_error_rate | {
  if [ "x$AWM_REPORT" == "xber" ]; then
    awk 'BEGIN { max_er = er = n = 0 } { er += $2; n++; if ($2 > max_er) max_er = $2;} END { print er / n, max_er; }'
  elif [ "x$AWM_REPORT" == "xfer" ]; then
    awk 'BEGIN { bad = n = 0 } { if ($2 > 0) bad++; n++; } END { print bad, n, bad * 100.0 / n; }'
  else
    echo "unknown report $AWM_REPORT" >&2
    exit 1
  fi
}
