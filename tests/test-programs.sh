#!/bin/bash

source test-common.sh

for TEST in testrawconverter
do
  if [ "x$Q" == "x1" ] && [ -z "$V" ]; then
    $TOP_BUILDDIR/src/$TEST > /dev/null
  else
    $TOP_BUILDDIR/src/$TEST
  fi
done

exit 0
