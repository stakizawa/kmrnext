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
! (type fccpx 1>/dev/null 2>/dev/null); fujitsu_env=$?
if test "$fujitsu_env" -eq "1"; then
   CC=fccpx
   CXX=FCCpx
fi

# Checks MPI and OpenMP
AC_LANG(C)
AX_MPI
AX_OPENMP([OPENMP_CFLAGS=$ax_cv_c_openmp])
AC_LANG(C++)
AX_MPI
AX_OPENMP([OPENMP_CXXFLAGS=$ax_cv_cxx_openmp])

# Checks other compiler flags
AC_SUBST(STD_CC)
AC_SUBST(STD_CXX)
AC_SUBST(WARN_CC)
AC_SUBST(WARN_CXX)

# The default CFLAGS and CXXFLAGS are remoted
case "$CC" in
gcc)
    CFLAGS=""
    CXXFLAGS=""
    STD_CC="-std=gnu99"
    STD_CXX="-std=c++98"
    WARN_CC="-Wall -Wextra -pedantic -Wno-long-long -Wno-unused-parameter -Wmissing-prototypes -Wshadow -Wconversion"
    WARN_CXX="-Wall -Wextra -pedantic -Wno-long-long -Wno-unused-parameter -Wshadow -Wconversion"
    ;;
fccpx)
    CFLAGS="-Kfast,preex -Xg"
    CXXFLAGS="-Kfast,preex -Xg"
    STD_CC="-std=gnu99"
    STD_CXX="-std=c++98"
    WARN_CC=""
    WARN_CXX=""
    ;;
esac

CFLAGS="-g $CFLAGS $WARN_CC $STD_CC $OPENMP_CFLAGS -MD"
CXXFLAGS="-g $CXXFLAGS $WARN_CXX $STD_CXX $OPENMP_CXXFLAGS -MD"

# Checks debug option
AC_ARG_ENABLE([debug],
	      [AS_HELP_STRING([--enable-debug],
			      [turn on debugging [default=no]])],
	      [], [enable_debug=no])
if test "$enable_debug" = "yes"; then
   CFLAGS="-DDEBUG -O0 $CFLAGS"
   CXXFLAGS="-DDEBUG -O0 $CXXFLAGS"
else
   CFLAGS="-O2 $CFLAGS"
   CXXFLAGS="-O2 $CXXFLAGS"
fi

# Checks KMR library [TODO]

# Checks KMR header [TODO]

# Checks for typedefs, structures, compiler characteristics,
# and library functions.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T
AC_FUNC_MALLOC

AC_OUTPUT
