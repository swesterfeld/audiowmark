#!/bin/bash

mkdir -p test
if [ ! -z "$(ls test)"  ]; then
  echo test dir not empty
  exit 1
fi

seq=1
cat test_list | while read f
do
  audiowmark gentest "$f" test/T$(printf "%02d__%s" $seq $(echo $f | sed 's, ,_,g;s,.*/,,g')).wav || exit 1
  ((seq++))
done
