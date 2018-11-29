#/bin/bash

echo 512 $(ber-test.sh)
for quality in 256 196 128 96 64
do
  echo $quality $(ber-test.sh double-mp3 $quality)
done
