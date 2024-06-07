#!/bin/bash
set -Eeo pipefail -x

# install dependencies
brew install autoconf-archive automake libsndfile fftw mpg123 libgcrypt libtool ffmpeg

# build zita-resampler
git clone https://github.com/swesterfeld/zita-resampler
cd zita-resampler
cmake .
sudo make install
cd ..
export DYLD_LIBRARY_PATH=/usr/local/lib:$DYLD_LIBRARY_PATH
# build audiowmark
./autogen.sh
make
make check
# build audiowmark with ffmpeg support
make clean
./autogen.sh --with-ffmpeg
make
make check
