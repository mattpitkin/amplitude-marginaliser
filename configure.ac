#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([amplitude-marginaliser], [0.1], [pitkin@gmail.com])

# add libtool
LT_INIT

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL
AC_PROG_INSTALL

AC_CONFIG_FILES([ \
  Makefile \
  src/Makefile \
  test/Makefile
])

# Checks for libraries.
AC_CHECK_LIB([m],[sin]) # math library

# check for gsl
AC_CHECK_LIB([gslcblas],[main])
AC_CHECK_LIB([gsl],[gsl_strerror],,[AC_MSG_ERROR([could not find the gsl library])])

# check for gsl headers
AC_CHECK_HEADERS([gsl/gsl_errno.h],,[AC_MSG_ERROR([could not find the gsl/gsl_errno.h header])])

# Checks for header files.
AC_CHECK_HEADERS([string.h sys/time.h])

# Checks for library functions.
AC_CHECK_FUNCS([gettimeofday memset])

# Initialise automake
AM_INIT_AUTOMAKE([1.10 no-define foreign])

AC_OUTPUT
