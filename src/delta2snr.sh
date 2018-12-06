for delta in 0.030 0.020 0.015 0.010 0.005 0.003 0.002 0.001
do
  #echo $delta $(AWM_PARAMS="--water-delta=$delta" ber-test.sh double-mp3 128) $(AWM_PARAMS="--water-delta=$delta" ber-test1.sh double-mp3 128)
  #echo $delta $(AWM_PARAMS="--water-delta=$delta" peaq.sh)
  echo $delta $(AWM_PARAMS="--water-delta=$delta" snr.sh)
done
