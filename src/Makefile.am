bin_PROGRAMS = audiowmark convcode

COMMON_SRC = utils.hh utils.cc

audiowmark_SOURCES = audiowmark.cc wavdata.cc wavdata.hh fft.cc fft.hh $(COMMON_SRC)
audiowmark_LDFLAGS = $(SNDFILE_LIBS) $(FFTW_LIBS)

convcode_SOURCES = convcode.cc $(COMMON_SRC)
