#pragma once

#include <acml.h>
#define BLAS(x,X) x
#define LAPACK(x,X) x
#define FC_GLOBAL(x,X) x

#define CblasColMajor
#define CblasTrans 'T'
#define CblasNoTrans 'N'
#define CblasRight 'R'
#define CblasLeft 'L'
#define CblasUpper 'U'

#define CLAPACK_INTERFACE

#define DGER(X1,X2,X3,X4,X5,X6,X7,X8,X9)                  BLAS(dger,DGER)(X1,X2,X3,X4,X5,X6,X7,X8,X9)
#define DSPMV(X1,X2,X3,X4,X5,X6,X7,X8,X9)                 BLAS(dspmv,DSPMV)(X1,X2,X3,X4,X5,X6,X7,X8,X9)
#define DTPMV(X1,X2,X3,X4,X5,X6,X7)                       BLAS(dtpmv,DTPMV)(X1,X2,X3,X4,X5,X6,X7)
#define DSYMM(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12)     BLAS(dsymm,DSYMM)(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12)
#define DGEMV(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11)         BLAS(dgemv,DGEMV)(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11)
#define DGEMM(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13) BLAS(dgemm,DGEMM)(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13)
#define DTRMM(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11)         BLAS(dtrmm,DTRMM)(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11)

#define DLANGE(X1,X2,X3,X4,X5,X6)       LAPACK(dlange,DLANGE)(X1,X2,X3,X4,X5)

#define DSPTRF(X1,X2,X3,X4,X5)          LAPACK(dsptrf,DSPTRF)(X1,X2,X3,X4,&X5)
#define DSPTRS(X1,X2,X3,X4,X5,X6,X7,X8) LAPACK(dsptrs,DSPTRS)(X1,X2,X3,X4,X5,X6,X7,&X8)
#define DSPTRI(X1,X2,X3,X4,X5,X6)       LAPACK(dsptri,DSPTRI)(X1,X2,X3,X4,&X6)
#define DPPTRF(X1,X2,X3,X4)             LAPACK(dpptrf,DPPTRF)(X1,X2,X3,&X4)
#define DPPTRI(X1,X2,X3,X4)             LAPACK(dpptri,DPPTRI)(X1,X2,X3,&X4)
#define DGETRF(X1,X2,X3,X4,X5)          { int Info; LAPACK(dgetrf,DGETRF)(X1,X2,X3,X4,X5,&Info); }
#define DGETRI(X1,X2,X3,X4)             { int Info; LAPACK(dgetri,DGETRI)(X1,X2,X3,X4,&Info);    }

#define DGESDD(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14) LAPACK(dgesdd,DGESDD)(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,&X14)
