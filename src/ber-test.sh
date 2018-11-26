# pseudo random pattern

PATTERN=4e1243bd22c66e76c2ba9eddc1f91394e57f9f83
TRANSFORM=$1

for i in test/T*
do
  echo $i
  audiowmark add $i t.wav $PATTERN >/dev/null
  if [ "x$TRANSFORM" == "xmp3" ]; then
    if [ "x$2" == "x" ]; then
      echo "need mp3 bitrate" >&2
      exit 1
    fi
    lame -b $2 t.wav t.mp3 --quiet
    rm t.wav
    ffmpeg -i t.mp3 t.wav -v quiet

    # some (low) mpeg quality settings use a lower sample rate
    if [ "x$(soxi -r t.wav)" != "x44100" ]; then
      sox t.wav tr.wav rate 44100
      mv tr.wav t.wav
    fi
  elif [ "x$TRANSFORM" == "x" ]; then
    :
  else
    echo "unknown transform $TRANSFORM" >&2
    exit 1
  fi
  audiowmark cmp $i t.wav $PATTERN
done | grep bit_error_rate | awk '{ er += $2; n++; if ($2 > max_er) max_er = $2;} END { print er / n, max_er; }'
