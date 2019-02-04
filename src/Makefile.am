bin_PROGRAMS = audiowmark

COMMON_SRC = utils.hh utils.cc convcode.hh convcode.cc random.hh random.cc mp3.cc wavdata.cc wavdata.hh
COMMON_LIBS = $(SNDFILE_LIBS) $(FFTW_LIBS) $(LIBGCRYPT_LIBS) $(LIBMPG123_LIBS)

audiowmark_SOURCES = audiowmark.cc fft.cc fft.hh $(COMMON_SRC)
audiowmark_LDFLAGS = $(COMMON_LIBS)

noinst_PROGRAMS = testconvcode testrandom testmp3

testconvcode_SOURCES = testconvcode.cc $(COMMON_SRC)
testconvcode_LDFLAGS = $(COMMON_LIBS)

testrandom_SOURCES = testrandom.cc $(COMMON_SRC)
testrandom_LDFLAGS = $(COMMON_LIBS)

testmp3_SOURCES = testmp3.cc $(COMMON_SRC)
testmp3_LDFLAGS = $(COMMON_LIBS)
