#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([KMR-Next], [1.0.0], [shinichiro.takizawa@riken.jp])
AC_CONFIG_SRCDIR([src/Makefile])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([config.make])

AC_CONFIG_MACRO_DIR([m4/])
m4_include([m4/ax_mpi.m4])
m4_include([m4/ax_openmp.m4])

# Checks compiler
AC_PROG_CC
AC_PROG_CXX
CFLAGS="$CFLAGS -Wall -MD"
CXXFLAGS="$CXXFLAGS -Wall -MD"

# Checks MPI and OpenMP
AC_LANG(C)
AX_MPI
AX_OPENMP([OPENMP_CFLAGS=$ax_cv_c_openmp])
AC_LANG(C++)
AX_MPI
AX_OPENMP([OPENMP_CXXFLAGS=$ax_cv_cxx_openmp])

CFLAGS="$CFLAGS $OPENMP_CFLAGS"
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"

# Checks debug option
AC_ARG_ENABLE([debug],
	      [AS_HELP_STRING([--enable-debug],
			      [turn on debugging [default=no]])],
	      [], [enable_debug=no])
if test "$enable_debug" = "yes"; then
   CFLAGS="-DDEBUG $CFLAGS"
   CXXFLAGS="-DDEBUG $CXXFLAGS"
fi


# Checks KMR library [TODO]

# Checks KMR header [TODO]

# Checks for typedefs, structures, compiler characteristics,
# and library functions.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_FUNC_MALLOC

AC_OUTPUT
