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

#include <sstream>
#include <cmath>
#include <cstdlib>

#include "MatLibConfig.h"
#include "matrix.h"
#include "symmatrix.h"
#include "om_utils.h"

SymMatrix::SymMatrix() : LinOp(0,0,SYMMETRIC,TWO),t(0),count(0) {}
SymMatrix::SymMatrix(const char* fname) : LinOp(0,0,SYMMETRIC,TWO),t(0),count(0) { this->load(fname); }
SymMatrix::SymMatrix(size_t N) : LinOp(N,N,SYMMETRIC,TWO),t(0),count(0) { alloc_data(); }
SymMatrix::SymMatrix(double* T, int* COUNT, size_t N) : LinOp(N,N,SYMMETRIC,TWO),t(T),count(COUNT) {(*count)++;}

bool SymMatrix::empty() const {return t==0;}
double* SymMatrix::data() const {return t;}
int* SymMatrix::DangerousGetCount() {return count;}

void SymMatrix::operator /=(double x) {(*this)*=(1/x);}

void SymMatrix::alloc_data()
{
    t=new double[(nlin()*(nlin()+1))/2];
    count=new int[1];
    (*count)=1;
}

void SymMatrix::destroy()
{
    if (t!=0) {
        (*count)--;
        if ((*count)==0) {
            delete[] t;
            delete[] count;
        }
    }
}

void SymMatrix::copy(const SymMatrix& A)
{
    t=A.t;
    nlin()=A.nlin();
    if (t) {
        count=A.count;
        (*count)++;
    }
}

SymMatrix SymMatrix::duplicate() const
{
    SymMatrix A(nlin());
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)(nlin()*(nlin()+1))/2,t,1,A.t,1);
#else
    for (size_t i=0;i<(nlin()*(nlin()+1))/2;i++)
        A.t[i]=t[i];
#endif
    return A;
}

const SymMatrix& SymMatrix::operator=(const SymMatrix& A) {
    destroy();
    copy(A);
    return *this;
}

const SymMatrix& SymMatrix::operator=(const double d) {
    for(size_t i=0;i<size();i++) t[i]=d;
    return *this;
}

SymMatrix::SymMatrix(const SymMatrix& A) {
    copy(A);
}

SymMatrix::SymMatrix(const Vector& v) {
    size_t N = v.size();
    nlin() = (size_t)((sqrt((double)(1+8*N))-1)/2+0.1);
    assert(nlin()*(nlin()+1)/2==N);
    t=v.data();
    if (t) {
        count=v.DangerousGetCount();
        (*count)++;
    }
}

SymMatrix::SymMatrix(const Matrix& A) {
    assert(A.nlin()==A.ncol());
    ncol() = A.ncol();
    alloc_data();
    for (size_t j=0; j<ncol(); j++)
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)(j+1),A.t+j*A.nlin(),1,t+(j*(j+1))/2,1);
#else
    for (size_t i=0; i<=j; i++)
        (*this)(i,j)=A(i,j);
#endif
}

void SymMatrix::set(double x) {
    for (size_t i=0;i<(nlin()*(nlin()+1))/2;i++)
        t[i]=x;
}

Vector SymMatrix::operator *(const Vector &v) const {
    assert(nlin()==v.size());
    Vector y(nlin());
#ifdef HAVE_BLAS
    DSPMV(CblasUpper,(int)nlin(),1.,t,v.t,1,0.,y.t,1);
#else
    for (size_t i=0;i<nlin();i++) {
        y(i)=0;
        for (size_t j=0;j<nlin();j++)
            y(i)+=(*this)(i,j)*v(j);
    }
#endif
    return y;
}

SymMatrix SymMatrix::inverse() const {
#ifdef HAVE_LAPACK
    SymMatrix invA=duplicate();
    // LU
    int *pivots=new int[nlin()];
    int info;
    DSPTRF('U',invA.nlin(),invA.t,pivots,info);
    // Inverse
    int size=(int)invA.nlin()*64;
    double *work=new double[size];
    DSPTRI('U',invA.nlin(),invA.t,pivots,work,info);

    delete[] pivots;
    delete[] work;
    return invA;
#else
    std::cerr << "!!!!! Inverse not implemented !!!!!" << std::endl;
    exit(1);
#endif
}

SymMatrix SymMatrix::posdefinverse() const {
    // supposes (*this) is definite positive
    SymMatrix invA=duplicate();
#ifdef HAVE_LAPACK
    // U'U factorization then inverse
    int info;
    DPPTRF('U',nlin(),invA.t,info);
    DPPTRI('U',nlin(),invA.t,info);
#else
    std::cerr << "Positive definite inverse not defined" << std::endl;
#endif
    return invA;
}

double SymMatrix::det() {
    SymMatrix invA=duplicate();
    double d = 1.0;
#ifdef HAVE_LAPACK
    // Bunch Kaufmqn
    int *pivots=new int[nlin()];
    int info;
    // TUDUtTt
    DSPTRF('U',invA.nlin(),invA.t,pivots,info);
    if(info<0)
        std::cout << "Big problem in det (DSPTRF)" << std::endl;
    for(int i = 0; i<(int) nlin(); i++){
        if(pivots[i] >= 0)
            d *= invA(i,i);
        else // pivots[i] < 0
            if(i < (int) nlin()-1 && pivots[i] == pivots[i+1]){
                d *= (invA(i,i)*invA(i+1,i+1)-invA(i,i+1)*invA(i+1,i));
                i++;
            }
            else
                std::cout << "Big problem in det" << std::endl;
    }
    delete[] pivots;
#else
    std::cerr << "Determinant not defined without LAPACK" << std::endl;
    exit(1);
#endif
    return(d);
}

void SymMatrix::eigen(Matrix & Z, Vector & D ){
    // performs the complete eigen-decomposition.
    //  (*this) = Z.D.Z'
    // -> eigenvector are columns of the Matrix Z.
    // (*this).Z[:,i] = D[i].Z[:,i]
#ifdef HAVE_LAPACK
    SymMatrix symtemp = duplicate();
    D = Vector(nlin());
    Z = Matrix(nlin(),nlin());

    int info;
    double lworkd;
    int lwork;
    int liwork;

    DSPEVD('V','U',nlin(),symtemp.t,D.t,Z.t,nlin(),&lworkd,-1,&liwork,-1,info);
    lwork = (int) lworkd;
    double * work = new double[lwork];
    int * iwork = new int[liwork];
    DSPEVD('V','U',nlin(),symtemp.t,D.t,Z.t,nlin(),work,lwork,iwork,liwork,info);

    delete[] work;
    delete[] iwork;
#endif
}

Matrix SymMatrix::operator *(const Matrix &B) const {
    assert(nlin()==B.nlin());
    Matrix C(nlin(),B.ncol());

#ifdef HAVE_BLAS
    Matrix D(*this);
    DSYMM(CblasLeft,  CblasUpper
        , (int)nlin(), (int)B.ncol(),
        1. , D.t, (int)D.ncol(),
        B.t, (int)B.nlin(),
        0, C.t,(int)C.nlin());
#else
    for (size_t j=0;j<B.ncol();j++)
        for (size_t i=0;i<ncol();i++)
        {
            C(i,j)=0;
            for (size_t k=0;k<ncol();k++)
                C(i,j)+=(*this)(i,k)*B(k,j);
        }
#endif
    return C;
}

SymMatrix SymMatrix::operator +(const SymMatrix &B) const {
    assert(nlin()==B.nlin());
    SymMatrix C=duplicate();
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(nlin()*(nlin()+1)/2), 1.0, B.t, 1, C.t , 1);
#else
    for (size_t i=0;i<nlin()*(nlin()+1)/2;i++)
        C.t[i]+=B.t[i];
#endif
    return C;
}

SymMatrix SymMatrix::operator -(const SymMatrix &B) const
{
    assert(nlin()==B.nlin());
    SymMatrix C=duplicate();
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(nlin()*(nlin()+1)/2), -1.0, B.t, 1, C.t , 1);
#else
    for (size_t i=0;i<nlin()*(nlin()+1)/2;i++)
        C.t[i]-=B.t[i];
#endif
    return C;
}

SymMatrix SymMatrix::operator *(double x) const {
    SymMatrix C(nlin());
    for (size_t k=0; k<nlin()*(nlin()+1)/2; k++) C.t[k] = t[k]*x;
    return C;
}

void SymMatrix::operator +=(const SymMatrix &B) {
    assert(nlin()==B.nlin());
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(nlin()*(nlin()+1)/2), 1.0, B.t, 1, t , 1);
#else
    for (size_t i=0;i<nlin()*(nlin()+1)/2;i++)
        t[i]+=B.t[i];
#endif
}

void SymMatrix::operator *=(double x) {
    for (size_t k=0; k<nlin()*(nlin()+1)/2; k++) t[k] *= x;
}

void SymMatrix::operator -=(const SymMatrix &B) {
    assert(nlin()==B.nlin());
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(nlin()*(nlin()+1)/2), -1.0, B.t, 1, t , 1);
#else
    for (size_t i=0;i<nlin()*(nlin()+1)/2;i++)
        t[i]+=B.t[i];
#endif
}

Matrix SymMatrix::operator()(size_t i_start, size_t i_end, size_t j_start, size_t j_end) const {
    Matrix retMat(i_end-i_start+1,j_end-j_start+1);
    for(size_t i=0;i<=i_end-i_start;i++)
        for(size_t j=0;j<=j_end-j_start;j++)
            retMat(i,j)=this->operator()(i_start+i,j_start+j);

    return retMat;
}

Matrix SymMatrix::submat(size_t istart, size_t isize, size_t jstart, size_t jsize) const {
    assert ( istart+isize<=nlin() && jstart+jsize<=nlin() );
    return (*this)(istart,istart+isize-1,jstart,jstart+jsize-1);
}

SymMatrix SymMatrix::submat(size_t istart, size_t iend) const {
    assert( iend > istart);
    size_t isize = iend - istart + 1;
    assert ( istart+isize<=nlin() );

    SymMatrix mat(isize);
    for(size_t i=istart;i<=iend;i++)
        for(size_t j=i;j<=iend;j++)
            mat(i,j)=this->operator()(i,j);

    return mat;
}

//returns the solution of (this)*X = B
Vector SymMatrix::solveLin(const Vector &B) const {
    SymMatrix invA=duplicate();
    Vector X = B.duplicate();

#ifdef HAVE_LAPACK
    // Bunch Kaufman Factorization
    int *pivots=new int[nlin()];
    int info;
    DSPTRF('U',invA.nlin(),invA.t,pivots,info);
    // Inverse
    int size=(int)invA.nlin()*64;
    double *work=new double[size];
    DSPTRS('U',invA.nlin(),1,invA.t,pivots,X.t,invA.nlin(),info);

    delete[] pivots;
    delete[] work;
#else
    std::cout << "solveLin not defined" << std::endl;
#endif
    return X;
}

// stores in B the solution of (this)*X = B, where B is a set of nbvect vector
void SymMatrix::solveLin(Vector * B, int nbvect) {
    SymMatrix invA=duplicate();

#ifdef HAVE_LAPACK
    // Bunch Kaufman Factorization
    int *pivots=new int[nlin()];
    int info;
    //char *uplo="U";
    DSPTRF('U',invA.nlin(),invA.t,pivots,info);
    // Inverse
    int size=(int) invA.nlin()*64;
    double *work=new double[size];
    for(int i = 0; i < nbvect; i++)
        DSPTRS('U',invA.nlin(),1,invA.t,pivots,B[i].t,invA.nlin(),info);

    delete[] pivots;
    delete[] work;
#else
    std::cout << "solveLin not defined" << std::endl;
#endif
}

void SymMatrix::info() const {
    if (nlin() == 0) {
        std::cout << "Matrix Empty" << std::endl;
        return;
    }

    std::cout << "Dimensions : " << nlin() << " x " << ncol() << std::endl;

    double minv = this->operator()(0,0);
    double maxv = this->operator()(0,0);
    size_t mini = 0;
    size_t maxi = 0;
    size_t minj = 0;
    size_t maxj = 0;

    for(size_t i = 0; i < nlin(); ++i)
    {
        for(size_t j = i; j < ncol(); ++j)
        {
            if (minv > this->operator()(i,j)) {
                minv = this->operator()(i,j);
                mini = i;
                minj = j;
            } else if (maxv < this->operator()(i,j)) {
                maxv = this->operator()(i,j);
                maxi = i;
                maxj = j;
            }
        }
    }
    std::cout << "Min Value : " << minv << " (" << mini << "," << minj << ")" << std::endl;
    std::cout << "Max Value : " << maxv << " (" << maxi << "," << maxj << ")" << std::endl;
    std::cout << "First Values" << std::endl;
    for(size_t i = 0; i < std::min(nlin(),(size_t) 5); ++i)
    {
        for(size_t j = i; j < std::min(ncol(),(size_t) 5); ++j)
        {
            std::cout << this->operator()(i,j) << " " ;
        }
        std::cout << std::endl ;
    }
}

// =======
// = IOs =
// =======

void SymMatrix::loadBin( const char *filename )
{
    maths::ifstream ifs(filename);
    ifs >> maths::format("binary") >> *this;
}

void SymMatrix::saveBin( const char *filename ) const
{
    maths::ofstream ofs(filename);
    ofs << maths::format("binary") << *this;
}

void SymMatrix::loadTxt( const char *filename )
{
    maths::ifstream ifs(filename);
    ifs >> maths::format("ascii") >> *this;
}

void SymMatrix::saveTxt( const char *filename ) const
{
    maths::ofstream ofs(filename);
    ofs << maths::format("ascii") << *this;
}

void SymMatrix::load( const char *filename ) {
    try {
        maths::ifstream ifs(filename);
        ifs >> *this;
    }
    catch (std::string s) {
        std::cout << s << std::endl;
    }
}

void SymMatrix::save( const char *filename ) const {
    try {
        maths::ofstream ofs(filename);
        ofs << *this;
    }
    catch (std::string s) {
        std::cout << s << std::endl;
    }
}

