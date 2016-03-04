#pragma once

#include <BlasLapackImplementations/FortranCInterface.h>

extern "C" {
    #include <cblas.h>
    #include <clapack.h>
}

#define BLAS(x,X) cblas_ ## x
#define LAPACK(x,X) clapack_ ## x

#define CLAPACK_INTERFACE

#include <BlasLapackImplementations/OM_C_BlasLapack.h>
#include <BlasLapackImplementations/OM_F_BlasLapack1.h>
