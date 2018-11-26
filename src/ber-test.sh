# pseudo random pattern

PATTERN=4e1243bd22c66e76c2ba9eddc1f91394e57f9f83

for i in test/T*
do
  echo $i
  audiowmark add $i t.wav $PATTERN >/dev/null
  audiowmark cmp $i t.wav $PATTERN
done | grep bit_error_rate | awk '{ er += $2; n++; if ($2 > max_er) max_er = $2;} END { print er / n, max_er; }'
