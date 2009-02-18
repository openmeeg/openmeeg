/* FILE: $Id$ */

/*
Project Name : OpenMEEG

version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Maureen.Clerc.AT.sophia.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.sophia.inria.fr)

The OpenMEEG software is a C++ package for solving the forward/inverse
problems of electroencephalography and magnetoencephalography.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's authors,  the holders of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

#ifndef OPENMEEG_MATLIBCONFIG_H
#define OPENMEEG_MATLIBCONFIG_H

//  cmake configuration.
#include <OpenMEEGConfigure.h>

#include "DLLDefinesOpenMEEGMaths.h"

#ifdef USE_MATIO
#include "matio.h"
#endif

#ifdef __APPLE__
#define F77_FUNC(name,NAME) name
#define F77_FUNC_(name,NAME) name ## _
#else
#define F77_FUNC(name,NAME) name ## _
#define F77_FUNC_(name,NAME) name ## _
#endif

#if defined(USE_ATLAS) || defined(USE_MKL)
#define HAVE_BLAS
#define HAVE_LAPACK
#endif

#ifdef HAVE_CONFIG_H
//  autoconf configuration.
#include <config.h>
#endif

#ifdef WIN32
    #pragma inline_recursion (on)
    #pragma inline_depth (255)
#endif

//#define inline __forceinline
//#define inline __attribute__((always_inline))
//#define inline __attribute__((weak)) inline

#if WIN32
    #pragma warning( disable : 4530)    //MSVC standard library can't be inlined
    #pragma warning( disable : 4996)    //MSVC warning C4996: declared deprecated
#endif

#ifdef USE_ATLAS
#ifdef __APPLE__
    extern "C" {
        #include <cblas.h>
        #include <clapack.h>
    }
    #define BLAS(x,X) cblas_ ## x
    #define LAPACK(x,X) x ## _
#else
    extern "C" {
        #include <atlas/cblas.h>
        #include <atlas/clapack.h>
    }
    #define BLAS(x,X) cblas_ ## x
    #define LAPACK(x,X) clapack_ ## x
#endif
#endif

#ifdef USE_MKL
// Hack to avoid the MKL declarations of Lapack Functions which do not use the power of C++ references
    #define _MKL_LAPACK_H_
    #include <mkl.h>
    #define BLAS(x,X) cblas_ ## x
    #define LAPACK(x,X) x
#endif

#ifdef USE_ACML
    #include <acml.h>
    //  Those macros are not used yet
    #define BLAS(x,X) x
    #define LAPACK(x,X) F77_FUNC(x,X)
    extern "C" void vrda_sin (int n, double *t, double *p);
    extern "C" void vrda_cos (int n, double *t, double *p);
    extern "C" void vrda_exp (int n, double *t, double *p);
    extern "C" void vrda_log (int n, double *t, double *p);
#endif

#if defined(HAVE_BLAS) && !defined(USE_ATLAS) && !defined(USE_ACML)
#ifndef USE_MKL
    #define CblasColMajor
    #define CblasTrans 'T'
    #define CblasNoTrans 'N'
    #define CblasRight 'R'
    #define CblasLeft 'L'
    #define CblasUpper 'U'
    #define BLAS(x,X) F77_FUNC(x,X)
    #define LAPACK(x,X) F77_FUNC(x,X)
#endif
    extern "C" {
#ifndef USE_MKL
        void BLAS(dcopy,DCOPY)(const int&,const double*,const int&,double*,const int);
        void BLAS(daxpy,DAXPY)(const int&,const double&,const double*,const int&,double*,const int&);
        double BLAS(ddot,DDOT)(const int&,const double*,const int&,const double*,const int&);
        double BLAS(dnrm2,DNRM2)(const int&,const double*,const int&);
        void BLAS(dscal,DSCAL)(const int&,const double&,double*,const int&);
        void BLAS(dspmv,DSPMV)(const char&,const int&,const double&,const double*,const double*,const int&,const double&,double*,const int&);
        void BLAS(dsymm,DSYMM)(const char&,const char&,const int&,const int&,const double&,const double*,const int&,const double*,const int&, const double*,double*,const int&);
        void BLAS(dgemm,DGEMM)(const char&,const char&,const int&,const int&,const int&,const double&,const double*,const int&,const double*,const int&,const double&,double*,const int&);
        void BLAS(dgemv,DGEMV)(const char&,const int&,const int&,const double&,const double*,const int&,const double*,const int&,const double&,double*,const int&);
#endif
        void LAPACK(dgetrf,DGETRF)(const int&,const int&,double*,const int&,int*,int&);
        void LAPACK(dgetri,DGETRI)(const int&,double*,const int&,int*,double*,const int&,int&);
    }
#endif

#if defined(HAVE_LAPACK)
    extern "C" {
        void F77_FUNC(dgesdd,DGESDD)(const char&,const int&,const int&,double*,const int&,double*,double*,const int&,double*,const int&,double*,const int&,int*,int&);
        void F77_FUNC(dpotf2,DPOTF2)(const char&,const int&,double*,const int&,int&);
        double F77_FUNC(dlange,DLANGE)(const char&,const int&,const int&,double*,const int&,double*);
        void F77_FUNC(dsptrf,DSPTRF)(const char&,const int&,double*,int*,int&);
        void F77_FUNC(dsptri,DSPTRI)(const char&,const int&,double*,int*,double*,int&);
        void F77_FUNC(dpptrf,DPPTRF)(const char&,const int&,double*,int&);
        void F77_FUNC(dpptri,DPPTRI)(const char&,const int&,double*,int&);
        void F77_FUNC(dspevd,DSPEVD)(const char&,const char&,const int&,double*,double*,double*,const int&,double*,const int&,int*,const int&,int&);
        void F77_FUNC(dsptrs,DSPTRS)(const char&,const int&,const int&,double*,int*,double*,const int&,int&);
    }
#endif

#define DGESDD F77_FUNC(dgesdd,DGESDD)
#define DPOTF2 F77_FUNC(dpotf2,DPOTF2)
#define DLANGE F77_FUNC(dlange,DLANGE)

#define DSPTRF F77_FUNC(dsptrf,DSPTRF)
#define DPPTRF F77_FUNC(dpptrf,DPPTRF)
#define DPPTRI F77_FUNC(dpptri,DPPTRI)
#define DSPEVD F77_FUNC(dspevd,DSPEVD)
#define DSPTRS F77_FUNC(dsptrs,DSPTRS)

#if defined(USE_ATLAS) || defined(USE_MKL)
    #define DSPMV(X1,X2,X3,X4,X5,X6,X7,X8,X9) BLAS(dspmv,DSPMV)(CblasColMajor,X1,X2,X3,X4,X5,X6,X7,X8,X9)
    #define DSYMM(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12) BLAS(dsymm,DSYMM)(CblasColMajor,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12)
    #define DGEMV(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11) BLAS(dgemv,DGEMV)(CblasColMajor,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11)
    #define DGEMM(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13) BLAS(dgemm,DGEMM)(CblasColMajor,X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13)
    #if defined(USE_ATLAS)
        #ifdef __APPLE__
            #define DGETRF(X1,X2,X3,X4,X5,X6) LAPACK(dgetrf,DGETRF)(&X1,&X2,X3,&X4,X5,&X6)
            // #define DGETRF(X1,X2,X3,X4,X5,X6) LAPACK(dgetrf,DGETRF)(CblasColMajor,X1,X2,X3,X4,X5)
            #define DGETRI(X1,X2,X3,X4,X5,X6,X7) LAPACK(dgetri,DGETRI)(&X1,X2,&X3,X4,X5,&X6,&X7)
            // #define DGETRI(X1,X2,X3,X4,X5,X6,X7) LAPACK(dgetri,DGETRI)(CblasColMajor,X1,X2,X3,X4)
        #else
            #define DGETRF(X1,X2,X3,X4,X5,X6) LAPACK(dgetrf,DGETRF)(CblasColMajor,X1,X2,X3,X4,X5)
            #define DGETRI(X1,X2,X3,X4,X5,X6,X7) LAPACK(dgetri,DGETRI)(CblasColMajor,X1,X2,X3,X4)
        #endif
    #else
        #define DGETRF(X1,X2,X3,X4,X5,X6) LAPACK(dgetrf,DGETRF)(X1,X2,X3,X4,X5,X6)
        #define DGETRI(X1,X2,X3,X4,X5,X6,X7) LAPACK(dgetri,DGETRI)(X1,X2,X3,X4,X5,X6,X7)
    #endif
    #define DSPTRI(X1,X2,X3,X4,X5,X6) F77_FUNC(dsptri,DSPTRI)(X1,X2,X3,X4,X5,X6)
#else
    #define DSPMV(X1,X2,X3,X4,X5,X6,X7,X8,X9) BLAS(dspmv,DSPMV)(X1,X2,X3,X4,X5,X6,X7,X8,X9)
    #define DSYMM(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12) BLAS(dsymm,DSYMM)(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12)
    #define DGEMV(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11) BLAS(dgemv,DGEMV)(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11)
    #define DGEMM(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13) BLAS(dgemm,DGEMM)(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13)
    #define DGETRF LAPACK(dgetrf,DGETRF)
    #if defined(USE_ACML)
        #define DGETRI(X1,X2,X3,X4,X5,X6,X7) LAPACK(dgetri,DGETRI)(X1,X2,X3,X4,X7)
        #define DSPTRI(X1,X2,X3,X4,X5,X6) LAPACK(dsptri,DSPTRI)(X1,X2,X3,X4,X6)
    #else
        #define DGETRI(X1,X2,X3,X4,X5,X6,X7) LAPACK(dgetri,DGETRI)(X1,X2,X3,X4,X5,X6,X7)
        #define DSPTRI(X1,X2,X3,X4,X5,X6) LAPACK(dsptri,DSPTRI)(X1,X2,X3,X4,X5,X6)
    #endif
#endif

#endif // !OPENMEEG_MATLIBCONFIG_H
