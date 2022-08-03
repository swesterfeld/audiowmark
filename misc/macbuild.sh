#!/bin/bash
set -Eeuo pipefail -x

brew install autoconf-archive automake libsndfile fftw mpg123 libgcrypt
./autogen.sh
make
make check
