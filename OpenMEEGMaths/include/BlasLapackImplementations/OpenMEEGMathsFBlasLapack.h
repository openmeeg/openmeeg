#pragma once

#if defined(USE_LAPACK)
    extern "C" {
        void LAPACK(dgesdd,DGESDD)(const char&,const int&,const int&,double*,const int&,double*,double*,const int&,double*,const int&,double*,const int&,int*,int&);
        void LAPACK(dpotf2,DPOTF2)(const char&,const int&,double*,const int&,int&);
        double LAPACK(dlange,DLANGE)(const char&,const int&,const int&,const double*,const int&,double*);
        void LAPACK(dsptrf,DSPTRF)(const char&,const int&,double*,int*,int&);
        void LAPACK(dtptri,DTPTRI)(const char&,const char&,const int&,double*,int&,int&);
        void LAPACK(dsptri,DSPTRI)(const char&,const int&,double*,int*,double*,int&);
        void LAPACK(dpptrf,DPPTRF)(const char&,const int&,double*,int&);
        void LAPACK(dpptri,DPPTRI)(const char&,const int&,double*,int&);
        void LAPACK(dspevd,DSPEVD)(const char&,const char&,const int&,double*,double*,double*,const int&,double*,const int&,int*,const int&,int&);
        void LAPACK(dsptrs,DSPTRS)(const char&,const int&,const int&,double*,int*,double*,const int&,int&);
    }
#endif

#define DGER  BLAS(dger,DGER)
#define DSPMV BLAS(dspmv,DSPMV)
#define DTPMV BLAS(dtpmv,DTPMV)
#define DSYMM BLAS(dsymm,DSYMM)
#define DGEMV BLAS(dgemv,DGEMV)
#define DGEMM BLAS(dgemm,DGEMM)
#define DTRMM BLAS(dtrmm,DTRMM)

#define DLANGE LAPACK(dlange,DLANGE)

#define DSPTRF LAPACK(dsptrf,DSPTRF)
#define DSPTRS LAPACK(dsptrs,DSPTRS)
#define DPPTRF LAPACK(dpptrf,DPPTRF)
#define DPPTRI LAPACK(dpptri,DPPTRI)

#define DGESDD FC_GLOBAL(dgesdd,DGESDD)

#define DGETRF LAPACK(dgetrf,DGETRF)

#define DGETRI LAPACK(dgetri,DGETRI)
#define DSPTRI LAPACK(dsptri,DSPTRI)
