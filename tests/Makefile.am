check: detect-speed-test block-decoder-test clip-decoder-test pipe-test short-payload-test

EXTRA_DIST = detect-speed-test.sh block-decoder-test.sh clip-decoder-test.sh pipe-test.sh short-payload-test.sh

detect-speed-test:
	Q=1 $(top_srcdir)/tests/detect-speed-test.sh

block-decoder-test:
	Q=1 $(top_srcdir)/tests/block-decoder-test.sh

clip-decoder-test:
	Q=1 $(top_srcdir)/tests/clip-decoder-test.sh

pipe-test:
	Q=1 $(top_srcdir)/tests/pipe-test.sh

short-payload-test:
	Q=1 $(top_srcdir)/tests/short-payload-test.sh
