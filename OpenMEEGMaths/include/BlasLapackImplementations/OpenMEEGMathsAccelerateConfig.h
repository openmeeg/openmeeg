#pragma once

// Accelerate ships no LAPACKE, so unlike OpenBLAS/MKL this goes through the
// Fortran interface. Its modern entry points need macOS >= 13.3; FC_ALIAS in
// FortranAlias.h is what actually selects them over the deprecated legacy ones.

#define FC_GLOBAL(x,X) x ## _

#define CblasColMajor
#define CblasTrans 'T'
#define CblasNoTrans 'N'
#define CblasRight 'R'
#define CblasLeft 'L'
#define CblasUpper 'U'

#define BLAS(x,X) FC_GLOBAL(x,X)
#define LAPACK(x,X) FC_GLOBAL(x,X)

// __LAPACK_int, in the LP64 interface we ask for (no ACCELERATE_LAPACK_ILP64)
typedef int BLAS_INT;

#define USE_LAPACK

#include <BlasLapackImplementations/blas.h>
#include <BlasLapackImplementations/lapack.h>
#include <BlasLapackImplementations/OpenMEEGMathsFBlasLapack.h>
