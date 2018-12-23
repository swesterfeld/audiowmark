SEEDS="$1"
MAX_SEED=$(($SEEDS - 1))
P1="$2"
P2="$3"
shift 3
echo "n seeds       : $SEEDS"
echo "ber-test args : $@"
echo "left          : $P1"
echo "right         : $P2"
for seed in $(seq 0 $MAX_SEED)
do
  echo $(AWM_SEEDS=$seed AWM_PARAMS="$P1" ber-test.sh "$@") $(AWM_SEEDS=$seed AWM_PARAMS="$P2" ber-test.sh "$@")
done | awk '{a += $1; if ($2 > b) b = $2; c += $3; if ($4 > d) d = $4; n++; } {printf ("%.5f %.5f     -     %.5f %.5f    -    (((%s)))\n", a/n, b, c/n, d, $0);}'
