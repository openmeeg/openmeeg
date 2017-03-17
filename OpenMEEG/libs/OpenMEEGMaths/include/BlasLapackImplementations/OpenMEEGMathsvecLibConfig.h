
#pragma once

#include <BlasLapackImplementations/FortranCInterface.h>

extern "C" {
    #include <cblas.h>
    #include <clapack.h>
}

#define BLAS(x,X) cblas_ ## x
#define LAPACK(x,X) x ## _

#define CLAPACK_INTERFACE

#include <BlasLapackImplementations/OpenMEEGMathsCBlasLapack.h>
#include <BlasLapackImplementations/OpenMEEGMathsFBlasLapack1.h>

#define DGETRF(X1,X2,X3,X4,X5,X6) LAPACK(dgetrf,DGETRF)(&X1,&X2,X3,&X4,X5,&X6)
#define DGETRI(X1,X2,X3,X4,X5,X6,X7) LAPACK(dgetri,DGETRI)(&X1,X2,&X3,X4,X5,&X6,&X7)
