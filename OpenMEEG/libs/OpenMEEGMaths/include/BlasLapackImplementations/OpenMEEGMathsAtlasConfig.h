#pragma once

#include <BlasLapackImplementations/FortranCInterface.h>

extern "C" {
    #include <cblas.h>
    #include <clapack.h>
}

#define BLAS(x,X) cblas_ ## x
#define LAPACK(x,X) clapack_ ## x

#define CLAPACK_INTERFACE

#include <BlasLapackImplementations/OpenMEEGMathsCBlasLapack.h>
#include <BlasLapackImplementations/OpenMEEGMathsFBlasLapack1.h>

#define DGETRF(X1,X2,X3,X4,X5) LAPACK(dgetrf,DGETRF)(CblasColMajor,X1,X2,X3,X4,X5)
#define DGETRI(X1,X2,X3,X4)    LAPACK(dgetri,DGETRI)(CblasColMajor,X1,X2,X3,X4)

