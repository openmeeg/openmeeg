#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT
AM_INIT_AUTOMAKE(OpenMEEG,0.1,Openmeeg-bugs@lists.gforge.inria.fr)
AC_CONFIG_SRCDIR([src/forward.cpp])
AM_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_INSTALL

F77=gfortran

AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_WRAPPERS

# Checks for libraries.

# BLAS, LAPACK, RPATH

#AC_F77_DUMMY_MAIN
AC_LANG(C)
ACX_BLAS()
#ACX_LAPACK

AC_CHECK_RPATH_OPTION  # !!! crash on mac os x !!!

if test acx_blas_ok; then
    ACX_LAPACK
    AC_SUBST(BLAS_INCLUDES)
    BLAS_INCLUDES=
    case $acx_blas_libs in
        *atlas*)
            BLAS_LIBS="${BLAS_LIBS} -lcblas"
            AC_DEFINE([USE_ATLAS],[],[Use the Atlas blas implementation])
            AC_CHECK_HEADER([cblas.h])
            AC_CHECK_HEADER([clapack.h])
            ;;
        -lmkl*)
            changequote(,)
            acx_blas_include=`echo $with_blas_dir | sed -e 's:/lib/[^/]*/\?$::'`/include
            changequote([,])
            AC_DEFINE([USE_MKL],[],[Use the mkl blas implementation])
            AC_CHECK_HEADER([$acx_blas_include/mkl.h])
            BLAS_INCLUDES="-I$acx_blas_include"
            BLAS_LIBS="$BLAS_LIBS @RPATH@$with_blas_dir"
            ;;
    esac
fi  

#if test acx_blas_ok; then
#    case $acx_blas_libs in
#        *atlas*)
#            BLAS_LIBS="${BLAS_LIBS} -lcblas"
#            AC_DEFINE([USE_ATLAS],[],[Use the Atlas blas implementation])
#            AC_CHECK_HEADER([cblas.h])
#            AC_CHECK_HEADER([clapack.h])
#            ;;
#    esac
#fi

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h limits.h stdlib.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt strrchr])

#                 Computations/Head1/Makefile
AC_CONFIG_FILES([
         Makefile
         src/Makefile
         libs/Makefile
         libs/MatLib/Makefile
         tests/Makefile
         libs/MatLib/Test/Makefile
])
AC_OUTPUT
