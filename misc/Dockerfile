FROM gcc:latest

RUN apt-get update && apt-get install -y \
  libgcrypt20-dev libsndfile1-dev libmpg123-dev libzita-resampler-dev \
  libfftw3-dev autoconf-archive clang

ADD . /audiowmark
WORKDIR /audiowmark

RUN misc/build.sh
