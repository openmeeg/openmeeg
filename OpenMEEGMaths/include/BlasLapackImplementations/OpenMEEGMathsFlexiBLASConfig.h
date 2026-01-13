#pragma once

#include <stdio.h>

#include <flexiblas/flexiblas_api.h>
#include <flexiblas/flexiblas_fortran_mangle.h>

#define CblasColMajor
#define CblasTrans 'T'
#define CblasNoTrans 'N'
#define CblasRight 'R'
#define CblasLeft 'L'
#define CblasUpper 'U'

#define BLAS(x,X) FC_GLOBAL(x,X)
#define LAPACK(x,X) FC_GLOBAL(x,X)

typedef FLEXIBLAS_API_INT BLAS_INT;

#define USE_LAPACK

#include <BlasLapackImplementations/blas.h>
#include <BlasLapackImplementations/lapack.h>
#include <BlasLapackImplementations/OpenMEEGMathsFBlasLapack.h>
