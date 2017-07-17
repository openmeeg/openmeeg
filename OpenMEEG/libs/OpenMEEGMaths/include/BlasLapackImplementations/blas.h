#pragma once

extern "C" {
    void BLAS(dcopy,DCOPY)(const int&,const double*,const int&,double*,const int&);
    void BLAS(daxpy,DAXPY)(const int&,const double&,const double*,const int&,double*,const int&);
    double BLAS(ddot,DDOT)(const int&,const double*,const int&,const double*,const int&);
    double BLAS(dnrm2,DNRM2)(const int&,const double*,const int&);
    void BLAS(dscal,DSCAL)(const int&,const double&,double*,const int&);
    void BLAS(dger,DGER)(const int&,const int&,const double&,const double*,const int&,const double*,const int&,double*,const int&);
    void BLAS(dspmv,DSPMV)(const char&,const int&,const double&,const double*,const double*,const int&,const double&,double*,const int&);
    void BLAS(dtpmv,DTPMV)(const char&,const char&,const char&,const int&,const double*,double*,const int&);
    void BLAS(dsymm,DSYMM)(const char&,const char&,const int&,const int&,const double&,const double*,const int&,const double*,const int&, const double&,double*,const int&);
    void BLAS(dgemm,DGEMM)(const char&,const char&,const int&,const int&,const int&,const double&,const double*,const int&,const double*,const int&,const double&,double*,const int&);
    void BLAS(dtrmm,DTRMM)(const char&,const char&,const char&,const char&,const int&,const int&,const double&,const double*,const int&,const double*,const int&);
    void BLAS(dgemv,DGEMV)(const char&,const int&,const int&,const double&,const double*,const int&,const double*,const int&,const double&,double*,const int&);
}
