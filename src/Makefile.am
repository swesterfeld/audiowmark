bin_PROGRAMS = audiowmark

COMMON_SRC = utils.hh utils.cc convcode.hh convcode.cc

audiowmark_SOURCES = audiowmark.cc wavdata.cc wavdata.hh fft.cc fft.hh $(COMMON_SRC)
audiowmark_LDFLAGS = $(SNDFILE_LIBS) $(FFTW_LIBS)

noinst_PROGRAMS = testconvcode

testconvcode_SOURCES = testconvcode.cc $(COMMON_SRC)
