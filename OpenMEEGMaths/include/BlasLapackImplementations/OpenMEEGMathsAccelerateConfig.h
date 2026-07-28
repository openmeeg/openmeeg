#pragma once

// Accelerate ships no LAPACKE, so unlike OpenBLAS/MKL this goes through the
// Fortran interface.

// Accelerate exports its modern (macOS >= 13.3) entry points as "dgemm$NEWLAPACK"
// rather than "dgemm_", so that suffix is simply the Fortran mangling we use;
// the plain names would silently bind Apple's frozen, deprecated LAPACK 3.2.1.
// The '$' is a compiler extension that -pedantic objects to, and it is only
// diagnosed here, where the token is lexed, not at each expansion.
#if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wdollar-in-identifier-extension"
#endif
#define FC_GLOBAL(x,X) x##$NEWLAPACK
#if defined(__clang__)
    #pragma clang diagnostic pop
#endif

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
