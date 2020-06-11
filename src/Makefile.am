bin_PROGRAMS = audiowmark
dist_bin_SCRIPTS = videowmark

COMMON_SRC = utils.hh utils.cc convcode.hh convcode.cc random.hh random.cc wavdata.cc wavdata.hh \
	     audiostream.cc audiostream.hh sfinputstream.cc sfinputstream.hh stdoutwavoutputstream.cc stdoutwavoutputstream.hh \
	     sfoutputstream.cc sfoutputstream.hh rawinputstream.cc rawinputstream.hh rawoutputstream.cc rawoutputstream.hh \
	     rawconverter.cc rawconverter.hh mp3inputstream.cc mp3inputstream.hh wmcommon.cc wmcommon.hh fft.cc fft.hh \
	     limiter.cc limiter.hh shortcode.cc shortcode.hh mpegts.cc mpegts.hh
COMMON_LIBS = $(SNDFILE_LIBS) $(FFTW_LIBS) $(LIBGCRYPT_LIBS) $(LIBMPG123_LIBS) $(FFMPEG_LIBS)

audiowmark_SOURCES = audiowmark.cc wmget.cc wmadd.cc $(COMMON_SRC)
audiowmark_LDFLAGS = $(COMMON_LIBS)

noinst_PROGRAMS = testconvcode testrandom testmp3 teststream testlimiter testshortcode testmpegts testhls

testconvcode_SOURCES = testconvcode.cc $(COMMON_SRC)
testconvcode_LDFLAGS = $(COMMON_LIBS)

testrandom_SOURCES = testrandom.cc $(COMMON_SRC)
testrandom_LDFLAGS = $(COMMON_LIBS)

testmp3_SOURCES = testmp3.cc $(COMMON_SRC)
testmp3_LDFLAGS = $(COMMON_LIBS)

teststream_SOURCES = teststream.cc $(COMMON_SRC)
teststream_LDFLAGS = $(COMMON_LIBS)

testlimiter_SOURCES = testlimiter.cc $(COMMON_SRC)
testlimiter_LDFLAGS = $(COMMON_LIBS)

testshortcode_SOURCES = testshortcode.cc $(COMMON_SRC)
testshortcode_LDFLAGS = $(COMMON_LIBS)

testmpegts_SOURCES = testmpegts.cc $(COMMON_SRC)
testmpegts_LDFLAGS = $(COMMON_LIBS)

testhls_SOURCES = testhls.cc wmadd.cc $(COMMON_SRC)
testhls_LDFLAGS = $(COMMON_LIBS)
