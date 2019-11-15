bin_PROGRAMS = audiowmark

COMMON_SRC = utils.hh utils.cc convcode.hh convcode.cc random.hh random.cc mp3.cc mp3.hh wavdata.cc wavdata.hh \
	     audiostream.cc audiostream.hh sfinputstream.cc sfinputstream.hh stdoutwavoutputstream.cc stdoutwavoutputstream.hh \
	     sfoutputstream.cc sfoutputstream.hh
COMMON_LIBS = $(SNDFILE_LIBS) $(FFTW_LIBS) $(LIBGCRYPT_LIBS) $(LIBMPG123_LIBS)

audiowmark_SOURCES = audiowmark.cc fft.cc fft.hh $(COMMON_SRC)
audiowmark_LDFLAGS = $(COMMON_LIBS)

noinst_PROGRAMS = testconvcode testrandom testmp3 teststream

testconvcode_SOURCES = testconvcode.cc $(COMMON_SRC)
testconvcode_LDFLAGS = $(COMMON_LIBS)

testrandom_SOURCES = testrandom.cc $(COMMON_SRC)
testrandom_LDFLAGS = $(COMMON_LIBS)

testmp3_SOURCES = testmp3.cc $(COMMON_SRC)
testmp3_LDFLAGS = $(COMMON_LIBS)

teststream_SOURCES = teststream.cc $(COMMON_SRC)
teststream_LDFLAGS = $(COMMON_LIBS)
