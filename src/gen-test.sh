#!/bin/bash

mkdir -p test

seq=1
cat test_list | while read f
do
  audiowmark gentest "$f" test/T$(printf "%02d__%s" $seq $(echo $f | sed 's, ,_,g;s,.*/,,g')).wav || exit 1
  ((seq++))
done
