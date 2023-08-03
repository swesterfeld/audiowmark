FROM gcc:latest

RUN apt-get update && apt-get install -y \
  build-essential automake autoconf libtool autoconf-archive gettext \
  libfftw3-dev libsndfile1-dev libgcrypt20-dev libzita-resampler-dev \
  libmpg123-dev

ADD . /audiowmark
WORKDIR /audiowmark

RUN ./autogen.sh
RUN make
RUN make install

VOLUME ["/data"]
WORKDIR /data

ENTRYPOINT ["/usr/local/bin/audiowmark"]
