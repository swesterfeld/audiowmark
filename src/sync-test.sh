#!/bin/bash
SEEDS="$1"
MAX_SEED=$(($SEEDS - 1))
P="$2"
shift 2
echo "n seeds       : $SEEDS"
echo "ber-test args : $@"
echo "params        : $P"
for seed in $(seq 0 $MAX_SEED)
do
  echo $(AWM_SEEDS=$seed AWM_PARAMS="$P" AWM_REPORT="sync" ber-test.sh "$@")
done | awk '{bad += $1; files += $2; print bad, files, bad * 100. / files }'
