# pseudo random pattern

PATTERN=4e1243bd22c66e76c2ba9eddc1f91394

for i in test/T*
do
  echo $i $(audiowmark add $i t.wav $PATTERN $AWM_PARAMS --snr 2>&1 | grep SNR)
done | awk '{s += $3; n++} END { print s/n; }'
