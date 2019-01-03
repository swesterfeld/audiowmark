bin_PROGRAMS = audiowmark

COMMON_SRC = utils.hh utils.cc convcode.hh convcode.cc random.hh random.cc

audiowmark_SOURCES = audiowmark.cc wavdata.cc wavdata.hh fft.cc fft.hh $(COMMON_SRC)
audiowmark_LDFLAGS = $(SNDFILE_LIBS) $(FFTW_LIBS) $(LIBGCRYPT_LIBS)

noinst_PROGRAMS = testconvcode testrandom

testconvcode_SOURCES = testconvcode.cc $(COMMON_SRC)
testconvcode_LDFLAGS = $(LIBGCRYPT_LIBS)

testrandom_SOURCES = testrandom.cc $(COMMON_SRC)
testrandom_LDFLAGS = $(LIBGCRYPT_LIBS)
