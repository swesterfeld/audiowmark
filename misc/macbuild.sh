#!/bin/bash
set -Eeuo pipefail -x

brew install autoconf-archive automake libsndfile fftw
./autogen.sh
make
make check
