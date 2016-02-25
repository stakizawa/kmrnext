#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([KMR-Next], [1.0.0], [shinichiro.takizawa@riken.jp])
AC_CONFIG_SRCDIR([src/Makefile])
AC_CONFIG_HEADERS([config.hpp])
AC_CONFIG_FILES([config.make])

AC_CONFIG_MACRO_DIR([m4/])
m4_include([m4/ax_mpi.m4])
m4_include([m4/ax_openmp.m4])

# Checks compiler
AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
! (type fccpx 1>/dev/null 2>/dev/null); fujitsu_env=$?
if test "$fujitsu_env" -eq "1"; then
    CC=fccpx
    CXX=FCCpx
    FC=frtpx
fi

# Checks MPI and OpenMP
AC_LANG(C)
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T
AX_MPI
AX_OPENMP([OPENMP_CFLAGS=$ax_cv_c_openmp])
AC_LANG(C++)
AC_TYPE_SIZE_T
AX_MPI
AX_OPENMP([OPENMP_CXXFLAGS=$ax_cv_cxx_openmp])
AC_LANG(Fortran)
AX_MPI
AX_OPENMP([OPENMP_FFLAGS=$ax_cv_fc_openmp])

# Checks other compiler flags
AC_SUBST(SYSENV)
AC_SUBST(STD_CC)
AC_SUBST(STD_CXX)
AC_SUBST(STD_FORT)
AC_SUBST(WARN_CC)
AC_SUBST(WARN_CXX)
AC_SUBST(WARN_FORT)
AC_SUBST(CXXLIB)

# The default CFLAGS, CXXFLAGS and FCFLAGS are removed
case "$CC" in
gcc)
    CFLAGS=""
    CXXFLAGS=""
    FCFLAGS="-cpp"
    STD_CC="-std=gnu99"
    STD_CXX="-std=c++98"
    STD_FORT="-std=gnu"
    WARN_CC="-Wall -Wextra -pedantic -Wno-long-long -Wno-unused-parameter -Wmissing-prototypes -Wshadow -Wconversion"
    WARN_CXX="-Wall -Wextra -pedantic -Wno-long-long -Wno-unused-parameter -Wshadow -Wconversion"
    WARN_FORT="-Wall -Wextra -pedantic -Wno-unused-parameter -Wintrinsic-shadow -Wconversion"
    CXXLIB="-lstdc++"
    _CXXLIB_MPI="-lmpi_cxx"
    SYSENV="GNULINUX"
    AC_DEFINE([SYSENV_GNULINUX], [1],
              [Define if system is gnu-linux])
    ;;
fccpx)
    CFLAGS="-Kfast,preex -Xg"
    CXXFLAGS="-Kfast,preex -Xg"
    FCFLAGS="-Kfast,preex"
    STD_CC="-std=gnu99"
    STD_CXX="-std=gnu++98"
    STD_FORT=""
    WARN_CC=""
    WARN_CXX=""
    WARN_FORT=""
    CXXLIB="-lfjdemgl -lstd_mt -lstdc++"
    _CXXLIB_MPI="-lmpi_cxx"
    SYSENV="K_FX"
    AC_DEFINE([SYSENV_K_FX], [1],
              [Define if system is K computer of FX10])
    ;;
esac

CFLAGS="-g $CFLAGS $WARN_CC $STD_CC $OPENMP_CFLAGS -MD"
CXXFLAGS="-g $CXXFLAGS $WARN_CXX $STD_CXX $OPENMP_CXXFLAGS -MD"
FCFLAGS="-g $FCFLAGS $WARN_FORT $STD_FORT $OPENMP_FFLAGS"
LDFLAGS="$OPENMP_CXXFLAGS"

# Checks debug option
AC_ARG_ENABLE([debug],
	      [AS_HELP_STRING([--enable-debug],
			      [turn on debugging [default=no]])],
	      [], [enable_debug=no])
if test "$enable_debug" = "yes"; then
    CFLAGS="-DDEBUG -O0 $CFLAGS"
    CXXFLAGS="-DDEBUG -O0 $CXXFLAGS"
    FCFLAGS="-DDEBUG -O0 $FCFLAGS"
else
    CFLAGS="-O2 $CFLAGS"
    CXXFLAGS="-O2 $CXXFLAGS"
    FCFLAGS="-O2 $FCFLAGS"
fi

# Checks the backend library
AC_ARG_WITH([backend],
	    [AS_HELP_STRING([--with-backend=@<:@serial/kmr@:>@],
			    [set the runtime backend [default=serial]])],
	    [with_backend=$withval],
	    [with_backend=serial])
if test "$with_backend" = "kmr"; then
    TARGET="KMR"
    AC_DEFINE([BACKEND_KMR], [1],
    	      [Define if kmr backend is set])
    CC=$MPICC
    CXX=$MPICXX
    FC=$MPIFC
    CXXLIB="$CXXLIB $_CXXLIB_MPI"
else
    # in case of "serial"
    TARGET="SERIAL"
fi
AC_SUBST(TARGET)

# Checks KMR library
AC_ARG_WITH([kmr],
	    [AS_HELP_STRING([--with-kmr=PATH],
			    [set path for kmr [default=/usr/local]])],
	    [KMRDIR=$withval],
	    [KMRDIR=/usr/local])
AC_SUBST(KMRDIR)

AC_OUTPUT
