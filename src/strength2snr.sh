for strength in 30 20 15 10 5 3 2 1
do
  echo $strength $(AWM_PARAMS="--strength=$strength" snr.sh)
done
