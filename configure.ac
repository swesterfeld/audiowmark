AC_INIT([audiowmark],[0.6.1])
AC_CONFIG_SRCDIR([src/audiowmark.cc])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([foreign])
AC_PROG_CXX
LT_INIT

dnl
dnl sndfile
dnl
AC_DEFUN([AC_SNDFILE_REQUIREMENTS],
[
    PKG_CHECK_MODULES(SNDFILE, [sndfile])
    AC_SUBST(SNDFILE_CFLAGS)
    AC_SUBST(SNDFILE_LIBS)
])

dnl
dnl libmpg123
dnl
AC_DEFUN([AC_LIBMPG123_REQUIREMENTS],
[
    PKG_CHECK_MODULES(LIBMPG123, [libmpg123])
    AC_SUBST(LIBMPG123_CFLAGS)
    AC_SUBST(LIBMPG123_LIBS)
])

dnl
dnl zita resampler
dnl
AC_DEFUN([AC_ZITA_REQUIREMENTS],
[
    AC_LIB_HAVE_LINKFLAGS([zita-resampler])
    if test "x$HAVE_LIBZITA_RESAMPLER" = "xno"; then
        AC_MSG_ERROR([You need to install libzita-resampler to build this package.])
    fi
])

dnl ffmpeg stuff
AC_DEFUN([AC_FFMPEG_REQUIREMENTS],
[
  PKG_CHECK_MODULES(FFMPEG, libavcodec libavformat libavutil libswresample)
])

dnl FFTW3
AC_DEFUN([AC_FFTW_CHECK],
[
    dnl this used to be optional, but is currently required
    PKG_CHECK_MODULES(FFTW, [fftw3f])

    SPECTMORPH_HAVE_FFTW=1
    if test $SPECTMORPH_HAVE_FFTW -gt 0; then
    fftw_save_CFLAGS="$CFLAGS"
    CFLAGS="$CFLAGS $FFTW_CFLAGS"
    AC_MSG_CHECKING([whether FFTW is recent enough])
      AC_COMPILE_IFELSE([
        AC_LANG_SOURCE([
          #include "fftw3.h"
          int x = FFTW_WISDOM_ONLY;
        ])
      ],[
        AC_MSG_RESULT([yes])
      ],[
        AC_MSG_RESULT([no])
        SPECTMORPH_HAVE_FFTW=0
      ])
    fi
    CFLAGS="$fftw_save_CFLAGS"
    AC_DEFINE_UNQUOTED(SPECTMORPH_HAVE_FFTW, $SPECTMORPH_HAVE_FFTW, [Whether libfftw3 is available])
])

AC_SNDFILE_REQUIREMENTS
AC_LIBMPG123_REQUIREMENTS
AC_ZITA_REQUIREMENTS
AC_FFTW_CHECK
AM_PATH_LIBGCRYPT

dnl -------------------- ffmpeg is optional ----------------------------
AC_ARG_WITH([ffmpeg], [AS_HELP_STRING([--with-ffmpeg], [build against ffmpeg libraries])], [], [with_ffmpeg=no])
if test "x$with_ffmpeg" != "xno"; then
  AC_FFMPEG_REQUIREMENTS
  HAVE_FFMPEG=1
else
  HAVE_FFMPEG=0
fi
AC_DEFINE_UNQUOTED(HAVE_FFMPEG, $HAVE_FFMPEG, [whether ffmpeg libs are available])
AM_CONDITIONAL([COND_WITH_FFMPEG], [test "x$with_ffmpeg" != "xno"])
dnl -------------------------------------------------------------------------

dnl -------------------- address sanitizer ----------------------------------
AC_ARG_ENABLE([asan], [AS_HELP_STRING([--enable-asan], [build using address sanitizer])],
[
  CXXFLAGS="$CXXFLAGS -fsanitize=address -fno-omit-frame-pointer -fstack-protector-all -fno-inline"
  LDFLAGS="$LDFLAGS -lasan"
],
[
  enable_asan=no
])
dnl -------------------------------------------------------------------------

dnl -------------------- undefined behaviour sanitizer ----------------------
AC_ARG_ENABLE([ubsan], [AS_HELP_STRING([--enable-ubsan], [build using undefined behaviour sanitizer])],
[
  CXXFLAGS="$CXXFLAGS -fsanitize=undefined -fno-omit-frame-pointer -fstack-protector-all -fno-inline"
  LDFLAGS="$LDFLAGS -lubsan"
],
[
  enable_ubsan=no
])
dnl -------------------------------------------------------------------------

dnl -------------------- thread sanitizer -----------------------------------
AC_ARG_ENABLE([tsan], [AS_HELP_STRING([--enable-tsan], [build using thread sanitizer])],
[
  CXXFLAGS="$CXXFLAGS -fsanitize=thread -fno-omit-frame-pointer -fstack-protector-all -fno-inline"
  LDFLAGS="$LDFLAGS -ltsan"
],
[
  enable_tsan=no
])
dnl -------------------------------------------------------------------------

dnl -------------------- leak sanitizer -------------------------------------
AC_ARG_ENABLE([lsan], [AS_HELP_STRING([--enable-lsan], [build using leak sanitizer])],
[
  CXXFLAGS="$CXXFLAGS -fsanitize=leak -fno-omit-frame-pointer -fstack-protector-all -fno-inline"
  LDFLAGS="$LDFLAGS -llsan"
],
[
  enable_lsan=no
])
dnl -------------------------------------------------------------------------

dnl -------------------- glibcxx assertions ----------------------------------
AC_ARG_ENABLE(debug-cxx,AS_HELP_STRING([--enable-debug-cxx], [setup compiler flags to do C++ STL debug checks]),
[
  CXXFLAGS="$CXXFLAGS -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC"
],
[
  enable_debug_cxx=no
])
dnl -------------------------------------------------------------------------

# need c++14 mode
AX_CXX_COMPILE_STDCXX_14(ext)

# use -Wall / -pthread if available
AC_LANG_PUSH([C++])
AX_CHECK_COMPILE_FLAG([-Wall],    [CXXFLAGS="$CXXFLAGS -Wall"])
AX_CHECK_COMPILE_FLAG([-pthread], [CXXFLAGS="$CXXFLAGS -pthread"])
AC_LANG_POP([C++])

# Less cluttered build output
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile tests/test-common.sh])
AC_OUTPUT

# Output summary message

echo
echo "---------------------------------------------------------------------------"
echo "$PACKAGE_NAME $PACKAGE_VERSION"
echo "---------------------------------------------------------------------------"
echo "  Use address sanitizer. . . . . . . .: $enable_asan (for debugging)"
echo "  Use undefined behaviour sanitizer. .: $enable_ubsan (for debugging)"
echo "  Use thread sanitizer . . . . . . . .: $enable_tsan (for debugging)"
echo "  Use leak sanitizer . . . . . . . . .: $enable_lsan (for debugging)"
echo "  Use stdc++ debug mode. . . . . . . .: $enable_debug_cxx (for debugging)"
echo "  Use ffmpeg libs. . . . . . . . . . .: $with_ffmpeg (required for HLS)"
