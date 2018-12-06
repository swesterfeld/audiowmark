# pseudo random pattern

PATTERN=4e1243bd22c66e76c2ba9eddc1f91394e57f9f83

for i in test/T*
do
  audiowmark add $i t.wav $PATTERN $AWM_PARAMS >/dev/null
  echo $i $(audiowmark snr $i t.wav)
done | grep snr | awk '{s += $3; n++} END { print s/n; }'
