#pragma once

#include <BlasLapackImplementations/FortranCInterface.h>

#define CblasColMajor
#define CblasTrans 'T'
#define CblasNoTrans 'N'
#define CblasRight 'R'
#define CblasLeft 'L'
#define CblasUpper 'U'

#define BLAS(x,X) FC_GLOBAL(x,X)
#define LAPACK(x,X) FC_GLOBAL(x,X)

#include <BlasLapackImplementations/blas.h>
#include <BlasLapackImplementations/lapack.h>
#include <BlasLapackImplementations/OpenMEEGMathsFBlasLapack.h>
