AUTOMAKE_OPTIONS = dist-zstd no-dist-gzip

# build & test everything (HLS included) during make distcheck
AM_DISTCHECK_CONFIGURE_FLAGS = --with-ffmpeg

SUBDIRS = src tests
ACLOCAL_AMFLAGS = -I m4

if COND_WITH_DOCS
  SUBDIRS += docs
endif

EXTRA_DIST = README.adoc Dockerfile
