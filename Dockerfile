FROM gcc:latest

RUN apt-get update && apt-get install -y build-essential
RUN apt-get install -y libfftw3-dev
RUN apt-get install -y libsndfile1-dev
RUN apt-get install -y automake
RUN apt-get install -y autoconf
RUN apt-get install -y libtool
RUN apt-get install -y autoconf-archive

ADD . /audiowmark
WORKDIR /audiowmark

RUN ./autogen.sh
RUN make
RUN make install

ENTRYPOINT ["/usr/local/bin/audiowmark"]
