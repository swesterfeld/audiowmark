CHECKS = detect-speed-test block-decoder-test clip-decoder-test \
       pipe-test short-payload-test sync-test sample-rate-test \
       key-test wav-pipe-test test-programs

if COND_WITH_FFMPEG
CHECKS += hls-test
endif

EXTRA_DIST = detect-speed-test.sh block-decoder-test.sh clip-decoder-test.sh \
       pipe-test.sh short-payload-test.sh sync-test.sh sample-rate-test.sh \
       key-test.sh hls-test.sh wav-pipe-test.sh test-programs.sh

check: $(CHECKS)

detect-speed-test:
	Q=1 $(top_srcdir)/tests/detect-speed-test.sh

block-decoder-test:
	Q=1 $(top_srcdir)/tests/block-decoder-test.sh

clip-decoder-test:
	Q=1 $(top_srcdir)/tests/clip-decoder-test.sh

pipe-test:
	Q=1 $(top_srcdir)/tests/pipe-test.sh

wav-pipe-test:
	Q=1 $(top_srcdir)/tests/wav-pipe-test.sh

short-payload-test:
	Q=1 $(top_srcdir)/tests/short-payload-test.sh

sync-test:
	Q=1 $(top_srcdir)/tests/sync-test.sh

sample-rate-test:
	Q=1 $(top_srcdir)/tests/sample-rate-test.sh

key-test:
	Q=1 $(top_srcdir)/tests/key-test.sh

hls-test:
	Q=1 $(top_srcdir)/tests/hls-test.sh

test-programs:
	Q=1 $(top_srcdir)/tests/test-programs.sh
