#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([KMR-Next], [1.0.0], [shinichiro.takizawa@riken.jp])
AC_CONFIG_SRCDIR([src/Makefile])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([config.make])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

CXXFLAGS="$CXXFLAGS -Wall -MD"
CFLAGS="$CFLAGS -Wall -MD"

# Checks for library functions.
AC_FUNC_MALLOC

AC_OUTPUT
