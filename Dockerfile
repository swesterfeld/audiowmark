FROM gcc:latest

RUN apt-get update && apt-get install -y build-essential
RUN apt-get install -y libfftw3-dev
RUN apt-get install -y libsndfile1-dev
RUN apt-get install -y automake
RUN apt-get install -y autoconf
RUN apt-get install -y libtool
RUN apt-get install -y autoconf-archive
RUN apt-get install -y libgcrypt20-dev
RUN apt-get install -y libzita-resampler-dev
RUN apt-get install -y libmpg123-dev

ADD . /audiowmark
WORKDIR /audiowmark

RUN ./autogen.sh
RUN make
RUN make install

VOLUME ["/data"]
WORKDIR /data

ENTRYPOINT ["/usr/local/bin/audiowmark"]
