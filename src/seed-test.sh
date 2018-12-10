MAX_SEED=$1
P1="$2"
P2="$3"
shift 3
echo "max seed      : $MAX_SEED"
echo "ber-test args : $@"
echo "left          : $P1"
echo "right         : $P2"
for seed in $(seq 0 $MAX_SEED)
do
  echo $(AWM_PARAMS="$P1 --seed $seed" ber-test.sh "$@") $(AWM_PARAMS="$P2 --seed $seed" ber-test.sh "$@")
done | awk '{a += $1; b += $2; c += $3; d += $4; n++; } {printf ("%.5f %.5f     -     %.5f %.5f    -    (((%s)))\n", a/n, b/n, c/n, d/n, $0);}'
