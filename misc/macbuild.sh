#!/bin/bash
set -Eeo pipefail -x

# install dependencies
brew install autoconf-archive automake libsndfile fftw mpg123 libgcrypt libtool ffmpeg@7
export PKG_CONFIG_PATH="$(brew --prefix ffmpeg@7)/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export PATH="$(brew --prefix ffmpeg@7)/bin:$PATH"

# build zita-resampler
git clone https://github.com/swesterfeld/zita-resampler
cd zita-resampler
cmake .
sudo make install
cd ..
export DYLD_LIBRARY_PATH=/usr/local/lib:$DYLD_LIBRARY_PATH

# build audiowmark
./autogen.sh
NPROC=`sysctl -n hw.ncpu`
make -j $NPROC
make  -j $NPROC check

# test build audiowmark with ffmpeg support
make clean
./autogen.sh --with-ffmpeg
make -j $NPROC

### unfortunately HLS is currently broken on macOS, so although it builds, make check will fail
###
###make -j $NPROC check
