#pragma once

#include <BlasLapackImplementations/FortranAlias.h>

extern "C" {
    void BLAS(dcopy,DCOPY)(const int&,const double*,const int&,double*,const int&) FC_ALIAS(dcopy);
    void BLAS(daxpy,DAXPY)(const int&,const double&,const double*,const int&,double*,const int&) FC_ALIAS(daxpy);
    double BLAS(ddot,DDOT)(const int&,const double*,const int&,const double*,const int&) FC_ALIAS(ddot);
    double BLAS(dnrm2,DNRM2)(const int&,const double*,const int&) FC_ALIAS(dnrm2);
    void BLAS(dscal,DSCAL)(const int&,const double&,double*,const int&) FC_ALIAS(dscal);
    void BLAS(dger,DGER)(const int&,const int&,const double&,const double*,const int&,const double*,const int&,double*,const int&) FC_ALIAS(dger);
    void BLAS(dspmv,DSPMV)(const char&,const int&,const double&,const double*,const double*,const int&,const double&,double*,const int&) FC_ALIAS(dspmv);
    void BLAS(dtpmv,DTPMV)(const char&,const char&,const char&,const int&,const double*,double*,const int&) FC_ALIAS(dtpmv);
    void BLAS(dsymm,DSYMM)(const char&,const char&,const int&,const int&,const double&,const double*,const int&,const double*,const int&,const double&,double*,const int&) FC_ALIAS(dsymm);
    void BLAS(dgemm,DGEMM)(const char&,const char&,const int&,const int&,const int&,const double&,const double*,const int&,const double*,const int&,const double&,double*,const int&) FC_ALIAS(dgemm);
    void BLAS(dtrmm,DTRMM)(const char&,const char&,const char&,const char&,const int&,const int&,const double&,const double*,const int&,const double*,const int&) FC_ALIAS(dtrmm);
    void BLAS(dgemv,DGEMV)(const char&,const int&,const int&,const double&,const double*,const int&,const double*,const int&,const double&,double*,const int&) FC_ALIAS(dgemv);
}
