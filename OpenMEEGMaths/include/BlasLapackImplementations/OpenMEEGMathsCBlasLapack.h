
#pragma once

#define DGER(X1,X2,X3,X4,X5,X6,X7,X8,X9)                    BLAS(dger,DGER)(CblasColMajor,X1,X2,X3,X4,X5,X6,X7,X8,X9)
#define DSPMV(X1,X2,X3,X4,X5,X6,X7,X8,X9)                   BLAS(dspmv,DSPMV)(CblasColMajor,X1,X2,X3,X4,X5,X6,X7,X8,X9)
#define DTPMV(X1,X2,X3,X4,X5,X6,X7)                         BLAS(dtpmv,DTPMV)(CblasColMajor,X1,X2,X3,X4,X5,X6,X7)
#define DSYMM(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12)       BLAS(dsymm,DSYMM)(CblasColMajor,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12)
#define DGEMV(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11)           BLAS(dgemv,DGEMV)(CblasColMajor,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11)
#define DGEMM(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13)   BLAS(dgemm,DGEMM)(CblasColMajor,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13)
#define DTRMM(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11)           BLAS(dtrmm,DTRMM)(CblasColMajor,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11)
