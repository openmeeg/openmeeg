#pragma once

// Hack to avoid the MKL declarations of Lapack Functions which do not use the power of C++ references
#define _MKL_LAPACK_H_

#include <mkl.h>
#define BLAS(x,X) cblas_ ## x
#define LAPACK(x,X) LAPACKE_ ## x

#define CLAPACK_INTERFACE

extern "C" {
    double dlange(const char&,const BLAS_INT&,const BLAS_INT&,double*,const int&,double*);
}

#define DLANGE dlange

#define DSPTRF(X1,X2,X3,X4,X5)          LAPACK(dsptrf,DSPTRF)(LAPACK_COL_MAJOR,X1,X2,X3,X4)
#define DSPTRS(X1,X2,X3,X4,X5,X6,X7,X8) LAPACK(dsptrs,DSPTRS)(LAPACK_COL_MAJOR,X1,X2,X3,X4,X5,X6,X7)
#define DSPTRI(X1,X2,X3,X4,X5,X6)       LAPACK(dsptri,DSPTRI)(LAPACK_COL_MAJOR,X1,X2,X3,X4)
#define DPPTRF(X1,X2,X3,X4)             LAPACK(dpptrf,DPPTRF)(LAPACK_COL_MAJOR,X1,X2,X3)
#define DPPTRI(X1,X2,X3,X4)             LAPACK(dpptri,DPPTRI)(LAPACK_COL_MAJOR,X1,X2,X3)
#define DGETRF(X1,X2,X3,X4,X5)          LAPACK(dgetrf,DGETRF)(LAPACK_COL_MAJOR,X1,X2,X3,X4,X5)
#define DGETRI(X1,X2,X3,X4)             LAPACK(dgetri,DGETRI)(LAPACK_COL_MAJOR,X1,X2,X3,X4)

#define DGESDD(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14) LAPACK(dgesdd,DGESDD)(LAPACK_COL_MAJOR,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10)

#include <BlasLapackImplementations/OpenMEEGMathsCBlasLapack.h>
