#pragma once

// Accelerate ships no LAPACKE, so unlike OpenBLAS/MKL this goes through the
// Fortran interface.

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

// Accelerate exports its modern (macOS >= 13.3) entry points as e.g.
// "_dgemm$NEWLAPACK", which is unspellable as a C identifier and cannot be an
// asm label either, since those must follow a declarator. Renaming here instead
// keeps blas.h/lapack.h shared with the other backends. Skipping this would
// silently bind Apple's frozen, deprecated LAPACK 3.2.1 symbols.
#if !defined(__PRAGMA_REDEFINE_EXTNAME)
    #error "Accelerate support needs #pragma redefine_extname (gcc/clang)"
#endif
// The '$' is what -pedantic objects to, and there is no spelling without it.
// This has to precede the #define, which is where the token gets lexed.
#if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wdollar-in-identifier-extension"
#endif

#define OM_STR_(s) #s
#define OM_STR(s) OM_STR_(s)
#define OM_NEWLAPACK(x) _Pragma(OM_STR(redefine_extname x##_ _##x##$NEWLAPACK))

OM_NEWLAPACK(dcopy)
OM_NEWLAPACK(daxpy)
OM_NEWLAPACK(ddot)
OM_NEWLAPACK(dnrm2)
OM_NEWLAPACK(dscal)
OM_NEWLAPACK(dger)
OM_NEWLAPACK(dspmv)
OM_NEWLAPACK(dtpmv)
OM_NEWLAPACK(dsymm)
OM_NEWLAPACK(dgemm)
OM_NEWLAPACK(dtrmm)
OM_NEWLAPACK(dgemv)

OM_NEWLAPACK(dgetrf)
OM_NEWLAPACK(dgetri)

OM_NEWLAPACK(dgesdd)
OM_NEWLAPACK(dpotf2)
OM_NEWLAPACK(dlange)
OM_NEWLAPACK(dsptrf)
OM_NEWLAPACK(dtptri)
OM_NEWLAPACK(dsptri)
OM_NEWLAPACK(dpptrf)
OM_NEWLAPACK(dpptri)
OM_NEWLAPACK(dspevd)
OM_NEWLAPACK(dsptrs)

#if defined(__clang__)
    #pragma clang diagnostic pop
#endif

#undef OM_NEWLAPACK
#undef OM_STR
#undef OM_STR_

#include <BlasLapackImplementations/blas.h>
#include <BlasLapackImplementations/lapack.h>
#include <BlasLapackImplementations/OpenMEEGMathsFBlasLapack.h>
