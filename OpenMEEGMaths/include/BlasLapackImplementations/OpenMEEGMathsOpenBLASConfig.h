#pragma once

// https://github.com/xianyi/OpenBLAS/wiki/How-to-use-OpenBLAS-in-Microsoft-Visual-Studio#visual-studio-2017-c2017-standard
// https://github.com/xianyi/OpenBLAS/issues/3661
// https://github.com/Reference-LAPACK/lapack/issues/683
// https://stackoverflow.com/questions/47520244/using-openblas-lapacke-in-visual-studio
#if defined(_MSC_VER)
    // https://developercommunity.visualstudio.com/t/_Fcomplex-and-_Dcomplex-removed-in-Windo/10108779#T-ND10301376
    #define _CRT_USE_C_COMPLEX_H
    #include "complex.h"
    #define LAPACK_COMPLEX_CUSTOM
    #define lapack_complex_float _Fcomplex
    #define lapack_complex_double _Dcomplex
#elif defined(__clang__)
    #pragma clang diagnostic ignored "-Wc99-extensions"
#endif

#include <cblas.h>
#include <lapacke.h>
#undef I // undefine this def due to complex.h that causes issues later

typedef int BLAS_INT;

#if defined(USE_SCIPY_OPENBLAS)
    #define BLAS(x,X) scipy_cblas_ ## x
    #define LAPACK(x,X) scipy_LAPACKE_ ## x
#else
    #define BLAS(x,X) cblas_ ## x
    #define LAPACK(x,X) LAPACKE_ ## x
#endif

#define CLAPACK_INTERFACE
#define UNUSED(expr) do{(void)(expr);}while(0)

#define DLANGE(X1,X2,X3,X4,X5,X6)       LAPACK(dlange,DLANGE)(LAPACK_COL_MAJOR,X1,X2,X3,X4,X5);UNUSED(X6)

#define DSPTRF(X1,X2,X3,X4,X5)          LAPACK(dsptrf,DSPTRF)(LAPACK_COL_MAJOR,X1,X2,X3,X4)
#define DSPTRS(X1,X2,X3,X4,X5,X6,X7,X8) LAPACK(dsptrs,DSPTRS)(LAPACK_COL_MAJOR,X1,X2,X3,X4,X5,X6,X7)
#define DSPTRI(X1,X2,X3,X4,X5,X6)       LAPACK(dsptri,DSPTRI)(LAPACK_COL_MAJOR,X1,X2,X3,X4)
#define DPPTRF(X1,X2,X3,X4)             LAPACK(dpptrf,DPPTRF)(LAPACK_COL_MAJOR,X1,X2,X3)
#define DPPTRI(X1,X2,X3,X4)             LAPACK(dpptri,DPPTRI)(LAPACK_COL_MAJOR,X1,X2,X3)
#define DGETRF(X1,X2,X3,X4,X5)          LAPACK(dgetrf,DGETRF)(LAPACK_COL_MAJOR,X1,X2,X3,X4,X5)
#define DGETRI(X1,X2,X3,X4)             LAPACK(dgetri,DGETRI)(LAPACK_COL_MAJOR,X1,X2,X3,X4)

#define DGESDD(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14) LAPACK(dgesdd,DGESDD)(LAPACK_COL_MAJOR,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10)

#include <BlasLapackImplementations/OpenMEEGMathsCBlasLapack.h>
