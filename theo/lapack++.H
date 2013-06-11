#ifndef LAPACKPP_H
#define LAPACKPP_H

#define F77_FUNC(x,X) x ## _

#include <iostream>
#include <algorithm>

namespace Maths {

    template <typename T>
    T max(const T t1,const T t2,const T t3) {
        return std::max(t1,std::max(t2,t3));
    }
}

extern "C" {

    //  Various forms of the SVD operation.

#   ifndef F77_FUNC
#   error "Unknown Fortran Compiler, check for a newer autoconf..."
#   endif

    //  Singular Value Decomposition.

#   define SGESVD77 F77_FUNC(sgesvd,SGESVD)
#   define DGESVD77 F77_FUNC(dgesvd,DGESVD)

    void SGESVD77(const char&,const char&,const int&,const int&,float*,const int&,float*,float*,const int&,float*,const int&,float*,const int&,int&);
    void DGESVD77(const char&,const char&,const int&,const int&,double*,const int&,double*,double*,const int&,double*,const int&,double*,const int&,int&);

    //  Generalized Singular Value Decomposition

#   define DGGSVD77 F77_FUNC(dggsvd,DGGSVD)

    void DGGSVD77(const char&,const char&,const char&,const int&,const int&,const int&,
                  int&,int&,double*,const int&,double*,const int&,double*,double*,
                  double*,const int&,double*,const int&,double*,const int&,double*,int*,int&);

#   define DGGLSE77 F77_FUNC(dgglse,DGGLSE)

    void DGGLSE77(const int&,const int&,const int&,double*,const int&,double*,const int&,double*,double*,double*,double*,const int&,int&);

#   define DGESVX77 F77_FUNC(dgesvx,DGESVX)

    void DGESVX77(const char&,const char&,const int&,const int&,double*,const int&,double*,
                  const int&,int*,const char&,double*,double*,double*,const int&,double*,
                  const int&,double&,double*,double*,double*,int*,int&);

#   define DGELS77 F77_FUNC(dgels,DGELS)

    void DGELS77(const char&,const int&,const int&,const int&,double*,const int&,double*,const int&,double*,const int&,int&);

    // Eigen elements of a symmetric matrix

#   define SSYEV77 F77_FUNC(ssyev,SSYEV)
#   define DSYEV77 F77_FUNC(dsyev,DSYEV)

    void SSYEV77(const char&,const char&,const int&,float*,const int&,float*,float*,const int&,int&);
    void DSYEV77(const char&,const char&,const int&,double*,const int&,double*,double*,const int&,int&);

    //  Eigen elements of a general matrix.

#   define SGEEV77 F77_FUNC(sgeev,SGEEV)
#   define DGEEV77 F77_FUNC(dgeev,DGEEV)

    void SGEEV77(const char&,const char&,const int&,float*,const int&,float*,float*,float*,const int&,float*,const int&,float*,const int&,int&);
    void DGEEV77(const char&,const char&,const int&,double*,const int&,double*,double*,double*,const int&,double*,const int&,double*,const int&,int&);

    //  Linear Least Square problem.
    //  Returns the minimum norm vector X that minimizes the L2 norm ||B-A*X||^2.

#   define SGELSS77 F77_FUNC(sgelss,SGELSS)
#   define DGELSS77 F77_FUNC(dgelss,DGELSS)

    void SGELSS77(const int&,const int&,const int&,float*,const int&,float*,const int&,float*,const float&,int&,float*,const int&,int&);
    void DGELSS77(const int&,const int&,const int&,double*,const int&,double*,const int&,double*,const double&,int&,double*,const int&,int&);

    //  LU decomposition.

#   define DGETRF77 F77_FUNC(dgetrf,DGETRF)
#   define DGETRI77 F77_FUNC(dgetri,DGETRI)

    void DGETRF77(const int&,const int&,double*,const int&,int*,int&);
    void DGETRI77(const int&,double*,const int&,const int*,double*,const int&,int&);

    //  QR RQ QL LQ decompositions.

#   define DGELQF77 F77_FUNC(dgelqf,DGELQF)
#   define DGEQLF77 F77_FUNC(dgeqlf,DGEQLF)
#   define DGEQRF77 F77_FUNC(dgeqrf,DGEQRF)
#   define DGERQF77 F77_FUNC(dgerqf,DGERQF)

    void DGELQF77(const int&,const int&,double*,const int&,double*,double*,const int&,int&);
    void DGEQLF77(const int&,const int&,double*,const int&,double*,double*,const int&,int&);
    void DGEQRF77(const int&,const int&,double*,const int&,double*,double*,const int&,int&);
    void DGERQF77(const int&,const int&,double*,const int&,double*,double*,const int&,int&);

#   define DGEQPF77 F77_FUNC(dgeqpf,DGEQPF)

    void DGEQPF77(const int&,const int&,double*,const int&,int*,double*,double*,int&);

#   define DORGLQ77 F77_FUNC(dorglq,DORGLQ)
#   define DORGQL77 F77_FUNC(dorgql,DORGQL)
#   define DORGQR77 F77_FUNC(dorgqr,DORGQR)
#   define DORGRQ77 F77_FUNC(dorgrq,DORGRQ)

    void DORGLQ77(const int&,const int&,const int&,double*,const int&,const double*,double*,const int&,int&);
    void DORGQL77(const int&,const int&,const int&,double*,const int&,const double*,double*,const int&,int&);
    void DORGQR77(const int&,const int&,const int&,double*,const int&,const double*,double*,const int&,int&);
    void DORGRQ77(const int&,const int&,const int&,double*,const int&,const double*,double*,const int&,int&);

    //  Cholesky Factorization of a Symmetric Positive Definite Matrix

#   define DPOTRF77 F77_FUNC(dpotrf,DPOTRF)

    void DPOTRF77(const char&,const int&,double*,const int&,int&); 

#   define DGESV77 F77_FUNC(dgesv,DGESV)

    void DGESV77(const int&,const int&,double*,const int&,const int*, 
                 double *,const int&,int&);
}

template <class T>
class Lapack {
};

//  Those routines are specific to C matrices... I would prefer using a real matrix class.
//  But for now this is the only alternative (TODO).

template <>
class Lapack<double> {
public:

    //  Solution of a system of n equations in n unknowns A*X=B.

    static void
    Solve(const unsigned n,double *A,double *B,double *X) {
        const int N = n;

        const int ldwork =3*n;

        double MF[n*n],dwork[ldwork],R[n],C[n],FErr[n],BErr[n],RCond;
        int    pivot[n],iwork[n];
        int    info;
        char   Equed;

        DGESVX77('E','N',n,1,A,n,MF,n,pivot,Equed,R,C,B,n,X,n,RCond,FErr,BErr,dwork,iwork,info);
    }

    //  Solution of a system of m equations in n unknowns A*X=C subject to the p constraints B*X=D.

    static void
    CstrLeastSqrs(const unsigned n,const unsigned m,const unsigned p,double *A,double *B,
                  double *C,double *D,double *X)
    {
        const int N = n;
        const int M = m;
        const int P = p;

        const int ldwork = N+P+Maths::max(M,N,P);

        double work[ldwork];
        int    info;

        DGGLSE77(N,M,P,A,M,B,P,C,D,X,work,ldwork,info);
    }
     
    static void
    LeastSqrs(const unsigned n,const unsigned m,const unsigned nrhs,double* A,double* B,
              double* Sing,const double RCond,unsigned& rank)
    {
        const int N    = n;
        const int M    = m;
        const int NRHS = nrhs;

        const int minmn = std::min(M,N);
        const int maxmn = std::max(M,N);

        const int ldwork = 3*minmn+Maths::max(2*minmn,NRHS,maxmn);

        double work[ldwork];
        int    info,Rank;

        DGELSS77(N,M,NRHS,A,N,B,maxmn,Sing,RCond,Rank,work,ldwork,info);

        rank = Rank;
    }

    //  The Svd routines have been arranged to appear to work directly on
    //  C matrices.

    static double*
    Svd(const unsigned n,const unsigned m,double* Mat,double *D) {
        const int N = n;
        const int M = m;

        const int minmn = std::min(m,n);
        const int maxmn = std::max(m,n);

        const int ldwork = std::max(3*minmn+maxmn,5*minmn);

        double *DUMMY=0,work[ldwork];
        int    info;

        std::cerr << M << ' ' <<N << ' ' << ldwork << ' ' << Mat << ' ' << D << std::endl;
        DGESVD77('O','N',M,N,Mat,M,D,DUMMY,1,DUMMY,1,work,ldwork,info);
        std::cerr <<"B" << std::endl;

        return D;
    }

    static void
    Svd(const unsigned n,const unsigned m,double* Mat) {
        double D[std::min(m,n)];
        Svd(n,m,Mat,D);
    }

    static double*
    SvdV(const unsigned n,const unsigned m,double* Mat,double *D,double *V) {
        const int N = n;
        const int M = m;
        const int minmn = std::min(M,N);
        const int maxmn = std::max(M,N);
        const int ldwork = std::max(3*minmn+maxmn,5*minmn);

        double *DUMMY=0,work[ldwork];
        int    info;

        DGESVD77('A','N',M,N,Mat,M,D,V,M,DUMMY,1,work,ldwork,info);

        return D;
    }

    static double*
    FullSvd(const unsigned n,const unsigned m,double* Mat,double *D,double *VT) {

        const int N = n;
        const int M = m;
        const int minmn = std::min(M,N);
        const int maxmn = std::max(M,N);
        const int ldwork = std::max(3*minmn+maxmn,5*minmn);

        double *DUMMY=0,work[ldwork];
        int    info;

        DGESVD77('S','O',M,N,Mat,M,D,VT,M,DUMMY,1,work,ldwork,info);

        return D;
    }

    static void
    SvdV(const unsigned n,const unsigned m,double* Mat,double *V) {
        double D[m];
        SvdV(n,m,Mat,D,V);
    }

    //  Eigenvalues of a symmetric matrix.

    static void
    SymmEigenVals(const unsigned n,double* mat,double* values) {
        int       info;
        const int ldwork = std::max(1U,3*n-1);
        double    work[ldwork];

        DSYEV77('N','U',n,mat,n,values,work,ldwork,info);
    }

    //  Eigenvectors of a symmetric matrix.
    //  On exit, mat is filled with the eigenvectors
    static void
    SymmEigenElmts(const unsigned n,double* mat,double* values) {
        int       info;
        const int ldwork = std::max(1U,3*n-1);
        double    work[ldwork];

        DSYEV77('V','U',n,mat,n,values,work,ldwork,info);
    }

    // Compute for an N-by-N real nonsymmetric matrix A, the eigenvalues and,
    // optionally, the left and/or right eigenvectors

    static void
    EigenVals(const unsigned n,double* mat,double* valuesReal,double* valuesIma) {
        double*   dummy;
        const int ldwork = std::max(1U,3*n);
        double    work[ldwork];
        int       info;

        DGEEV77('N','N',n,mat,n,valuesReal,valuesIma,dummy,1,dummy,1,work,ldwork,info);
    }

    //  Eigenvectors of a non symmetric matrix.

    static void
    EigenElmts(const unsigned n,double* mat,double* valuesReal,double* valuesIma,double* vects) {
        double*   dummy;
        const int ldwork = 4*n;
        double    work[ldwork];
        int       info;
        DGEEV77('N','V',n,mat,n,valuesReal,valuesIma,dummy,1,vects,n,work,ldwork,info);
    }

    //  Generalized SVD 
    //  U^t * A * Q = SIGa   dimA=[MxN]
    //  V^t * B * Q = SIGb   dimB=[PxN]
    //  svd computes only Q [N*N] and SIGa [N] SIGb[N]
    //  CAUTION: A and B are modified: make a copy if you want to keep them

    static void
    gsvd(const unsigned m,const unsigned n,const unsigned p, 
         double *A,double* B,double *Q,double* U,double* V,double* SigA,double* SigB,int & K,int & L) {
        
        double *DUMMY=0;
        const int M = m,P = p,N = n; 
        const unsigned ws = Maths::max(M,P,3*N)+N;
        double work1[ws];
        int iwork[N];
        int info;
        DGGSVD77('U','V','Q',M,N,P,K,L,A,M,B,P,SigA,SigB,U,M,V,P,Q,N,work1,iwork,info);
    }

    static void
    gsvd(const unsigned m,const unsigned n,const unsigned p, 
         double *A,double* B,double *Q,double* SigA,double* SigB,int & K,int & L) {
        
        double *DUMMY=0;
        const int M = m,P = p,N = n; 
        const unsigned ws = Maths::max(M,P,3*N)+N;
        double work1[ws];
        int iwork[N];
        int info;
        DGGSVD77('N','N','Q',M,N,P,K,L,A,M,B,P,SigA,SigB,DUMMY,1,DUMMY,1,Q,N,work1,iwork,info);
    }

    //  Cholesky Decomposition of a Symmetric Positive Definite Matrix
    //  Decompose A into M^T M such as M is triangular
    //  Input: Mat is the whole symmetric matrix of size [nxn]    
    //  Output: Mat contains the upper-triangular matrix M

    static bool 
    cholesky(const unsigned n,double *Mat) {

        for (int i=0; i<n; i++) 
            for (int j=i+1; j<n; j++) Mat[i*n+j]=0.0;
        
        int info;
        DPOTRF77('L',n,Mat,n,info); 
        return (info==0);
    }

    static void
    LU(const unsigned m,const unsigned n,double *Mat,int* pivot) {
        int info;
        DGETRF77(m,n,Mat,m,pivot,info);
    }

    static double
    det(const unsigned n,double *Mat) {
        // BUG in FC5 x86_84 environment and atlas libs...
        // Apparently info is not reset to 0 even though the sources state that (a compiler bug) ??
        // Long term should we call directly the atlas routines ??
        
        int info = 0;
        int pivot[n];

        DGETRF77(n,n,Mat,n,pivot,info);

        if (info>0)
            return 0.0;
        
        double det = 1.0;
        for (int i=0;i<n;i++)
            det *= (pivot[i]!=(i+1)?-1:1);
	
        for (int i=0,l=0;i<n;++i,l+=n+1)
            det *= Mat[l];
        
        return det;
    }

    static void
    inverse(const unsigned n,double *Mat) {
        int    info;
        int    pivot[n];
        double work[n];

        DGETRF77(n,n,Mat,n,pivot,info);
        DGETRI77(n,Mat,n,pivot,work,n,info);
    }
 
    //  All these are working on fortran matrices. Correct this (TODO).

    static void
    QL(const unsigned m,const unsigned n,double *Mat,double *L) {
        const int M = m;
        const int N = n;
        const int ldwork = std::max(1,N);
        const int minmn  = std::min(M,N);

        double work[ldwork],Tau[minmn];
        int    info;

        DGEQLF77(M,N,Mat,M,Tau,work,ldwork,info);

        for (int i=0;i<N;++i)
            for (int j=i,k=i+(M-N);j<N;++j,++k)
                L[j*N+i] = Mat[i*M+k];

        DORGQL77(M,N,minmn,Mat,M,Tau,work,ldwork,info);
    }

    static void
    LQ(const unsigned m,const unsigned n,double *Mat,double *L) {
        const int M = m;
        const int N = n;
        const int ldwork = std::max(1,M);
        const int minmn  = std::min(M,N);

        double work[ldwork],Tau[minmn];
        int    info;

        DGELQF77(M,N,Mat,M,Tau,work,ldwork,info);

        for (int  i=0;i<M;++i)
            for (int  j=0;j<std::min(N,i+1);++j)
                L[i*M+j] = Mat[j*M+i];

        DORGLQ77(M,N,minmn,Mat,M,Tau,work,ldwork,info);
    }

    static void
    RQ(const unsigned m,const unsigned n,double *Mat,double *R) {
        const int M = m;
        const int N = n;
        const int ldwork = std::max(1,M);

        double work[ldwork],Tau[std::min(M,N)];
        int    info;

        DGERQF77(M,N,Mat,M,Tau,work,ldwork,info);

        for (int j=0;j<N;++j)
            for (int i=0;i<=j+(M-N);++i)
                R[j*M+i] = Mat[j*M+i];

        DORGRQ77(N,N,N,Mat+M-N,M,Tau,work,ldwork,info);
    }

    static void
    QR(const unsigned m,const unsigned n,double *Mat,double *R) {
        const int M = m;
        const int N = n;
        const int ldwork = Maths::max(1,M,N);
        const int minmn = std::min(M,N);

        double work[ldwork],Tau[minmn];
        int    info;

        DGEQRF77(M,N,Mat,M,Tau,work,ldwork,info);

        for (int i=0;i<N;++i)
            for (int j=0;j<=i;++j)
                R[j*N+i] = Mat[i*M+j];

        DORGQR77(M,N,minmn,Mat,M,Tau,work,ldwork,info);
    }
};

#endif // !LAPACKPP_H
