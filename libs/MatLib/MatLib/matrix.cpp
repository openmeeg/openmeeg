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

#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <cfloat>

#include "om_utils.h"
#include "matrix.h"
#include "symmatrix.h"
#include "sparse_matrix.h"
#include "vector.h"

const Matrix& Matrix::set(const double d) {
    for(size_t i=0;i<size();i++) data()[i]=d;
    return *this;
}

Matrix::Matrix(const SymMatrix& A): LinOp(A.nlin(),A.ncol(),FULL,TWO),value(new LinOpValue(size())) {
    for (size_t j=0; j<ncol();++j)
        for (size_t i=0; i<nlin();++i)
            (*this)(i,j) = A(i,j);
}

Matrix::Matrix(const Vector& v,const size_t M,const size_t N): LinOp(M,N,FULL,TWO) {
    assert(M*N==v.size());
    value = v.value;
}

Vector Matrix::operator*(const Vector &v) const
{
    assert(ncol()==v.nlin());
    Vector y(nlin());
#ifdef HAVE_BLAS
    DGEMV(CblasNoTrans,(int)nlin(),(int)ncol(),1.0,data(),(int)nlin(),v.data(),1,0.,y.data(),1);
#else
    for (size_t i=0;i<nlin();i++) {
        y(i) = 0;
        for (size_t j=0;j<ncol();j++)
            y(i) += (*this)(i,j)*v(j);
    }
#endif

    return y;
}

Matrix Matrix::submat(size_t istart, size_t isize, size_t jstart, size_t jsize) const {
    assert (istart+isize<=nlin() && jstart+jsize<=ncol());

    Matrix a(isize,jsize);
    for (size_t j=0; j<jsize; j++)
#ifdef HAVE_BLAS
        BLAS(dcopy,DCOPY)((int)(isize),data()+istart+(jstart+j)*nlin(),1,a.data()+j*isize,1);
#elif USE_ACML
        dcopy((int)(isize),data()+istart+(jstart+j)*nlin(),1,a.data()+j*isize,1);
#else
        for (size_t i=0; i<isize; i++)
            a(i,j) = (*this)(istart+i,jstart+j);
#endif
    return a;
}

Vector Matrix::getcol(size_t j) const {
    assert(j<ncol());
    Vector v(nlin());
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)nlin(),data()+nlin()*j,1,v.data(),1);
#else
    for (size_t i=0;i<nlin();i++) v.data()[i]=data()[i+nlin()*j];
#endif
    return v;
}

Vector Matrix::getlin(size_t i) const {
    assert(i<nlin());
    Vector v(ncol());
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)ncol(),data()+i,(int)nlin(),v.data(),1);
#else
    for (size_t j=0;j<ncol();j++) v.data()[j]=data()[i+nlin()*j];
#endif
    return v;
}

void Matrix::setcol(size_t j,const Vector& v) {
    assert(v.size()==nlin() && j<ncol());
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)nlin(),v.data(),1,data()+nlin()*j,1);
#else
    for (size_t i=0;i<nlin();i++) data()[i+nlin()*j]=v.data()[i];
#endif
}

void Matrix::setlin(size_t i,const Vector& v) {
    assert(v.size()==ncol() && i<nlin());
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)ncol(),v.data(),1,data()+i,(int)nlin());
#else
    for (size_t j=0;j<ncol();j++) data()[i+nlin()*j]=v.data()[j];
#endif
}

Vector Matrix::tmult(const Vector &v) const
{
    assert(nlin()==v.nlin());
    Vector y(ncol());
#ifdef HAVE_BLAS
    DGEMV(CblasTrans,(int)nlin(),(int)ncol(),1.,data(),(int)nlin(),v.data(),1,0.,y.data(),1);
#else
    for (size_t i=0;i<ncol();i++) {
        y(i)=0;
        for (size_t j=0;j<nlin();j++)
            y(i)+=(*this)(j,i)*v(j);
    }
#endif

    return y;
}

Matrix Matrix::inverse() const
{
#ifdef HAVE_LAPACK
    assert(nlin()==ncol());
    Matrix invA(*this,DEEP_COPY);
    // LU
    int *pivots=new int[ncol()];
    int info;
    int nlin_local = invA.nlin();
    int nlin_local2 = invA.nlin();
    int ncol_local = invA.ncol();
    DGETRF(nlin_local,ncol_local,invA.data(),nlin_local2,pivots,info);
    // DGETRF(invA.nlin(),invA.ncol(),invA.data(),invA.nlin(),pivots,info);
    // Inverse
    int size=(int)invA.ncol()*64;
    double *work=new double[size];
    DGETRI(ncol_local,invA.data(),ncol_local,pivots,work,size,info);
    // DGETRI(invA.ncol(),invA.data(),invA.ncol(),pivots,work,size,info);
    delete[] pivots;
    delete[] work;
    return invA;
#else
    std::cerr << "!!!!! Inverse not implemented !!!!!" << std::endl;
    exit(1);
#endif
}

Matrix Matrix::pinverse(double tolrel) const {
#if defined(HAVE_BLAS) && defined(HAVE_LAPACK)
    if(ncol() > nlin()) return transpose().pinverse().transpose();
    else {
        Matrix result(ncol(),nlin());
        Matrix U,S,V;
        svd(U,S,V);
        double maxs=0;
        int mimi=(int)std::min(S.nlin(),S.ncol());
        for(int i=0;i<mimi;i++) maxs=std::max(S(i,i),maxs);
        if (tolrel==0) tolrel=DBL_EPSILON;
        double tol = std::max(nlin(),ncol()) * maxs * tolrel;
        int r=0; for(int i=0;i<mimi;i++) if(S(i,i)>tol) r++;
        if (r == 0) {
            result.set(0.);
            return result;
        } else {
            Matrix s(r,r); s.set(0);
            for(int i=0;i<r;i++) s(i,i)=1.0/S(i,i);
            const Matrix Vbis(V,r);
            const Matrix Ubis(U,r);
            return Vbis*s*Ubis.transpose();
        }
    }
#else
    std::cerr << "pinv not implemented without blas/lapack" << std::endl;
    exit(1);
#endif
}

Matrix Matrix::transpose() const {
    Matrix result(ncol(),nlin());
    for(size_t i=0;i<nlin();i++) for(size_t j=0;j<ncol();j++) result(j,i)=(*this)(i,j);
    return result;
}

void Matrix::svd(Matrix &U,Matrix &S, Matrix &V) const {
#ifdef HAVE_LAPACK
    Matrix cpy(*this,DEEP_COPY);
    int mimi = (int)std::min(nlin(),ncol());
    U = Matrix(nlin(),ncol()); U.set(0);
    V = Matrix(ncol(),ncol()); V.set(0);
    S = Matrix(ncol(),ncol()); S.set(0);
    double *s=new double[mimi];
    int lwork=4 *mimi*mimi + (int)std::max(nlin(),ncol()) + 9*mimi;
    double *work=new double[lwork];
    int *iwork=new int[8*mimi];
    int info;
    DGESDD('S',nlin(),ncol(),cpy.data(),nlin(),s,U.data(),U.nlin(),V.data(),V.nlin(),work,lwork,iwork,info);
    for(int i=0;i<mimi;i++) S(i,i)=s[i];
    V=V.transpose();
    delete[] s;
    delete[] work;
    delete[] iwork;
#else
    std::cerr<<"svd not implemented without blas/lapack"<<std::endl;
#endif
}

Matrix Matrix::operator *(const Matrix &B) const
{
    assert(ncol()==B.nlin());
    size_t p=ncol();
    Matrix C(nlin(),B.ncol());
#ifdef HAVE_BLAS
    DGEMM(CblasNoTrans,CblasNoTrans,
        (int)C.nlin(),(int)C.ncol(),(int)p,
        1.,data(),(int)nlin(),
        B.data(),(int)B.nlin(),
        0.,C.data(),(int)C.nlin());
#else
    for (size_t i=0;i<C.nlin();i++)
        for (size_t j=0;j<C.ncol();j++) {
            C(i,j)=0;
            for (size_t k=0;k<p;k++)
                C(i,j)+=(*this)(i,k)*B(k,j);
        }
#endif
        return C;
}

Matrix Matrix::operator *(const SparseMatrix &mat) const
{
    assert(ncol()==mat.nlin());
    Matrix out(nlin(),mat.ncol());
    out.set(0.0);

    SparseMatrix::const_iterator it;
    for(it = mat.begin(); it != mat.end(); ++it) {
        size_t i = it->first.first;
        size_t j = it->first.second;
        double val = it->second;
        for(size_t k = 0; k < nlin(); ++k) {
            out(k,j) += this->operator()(k,i) * val;
        }
    }
    return out;
}

Matrix Matrix::tmult(const Matrix &B) const
{
    assert(nlin()==B.nlin());
    size_t p=nlin();
    Matrix C(ncol(),B.ncol());
#ifdef HAVE_BLAS
    DGEMM(CblasTrans,CblasNoTrans,
        (int)C.nlin(),(int)C.ncol(),(int)p,
        1.,data(),(int)nlin(),
        B.data(),(int)B.nlin(),
        0.,C.data(),(int)C.nlin());
#else
    for (size_t i=0;i<C.nlin();i++)
        for (size_t j=0;j<C.ncol();j++) {
            C(i,j)=0;
            for (size_t k=0;k<p;k++)
                C(i,j)+=(*this)(k,i)*B(k,j);
        }
#endif
        return C;
}

Matrix Matrix::multt(const Matrix &B) const
{
    assert(ncol()==B.ncol());
    size_t p=ncol();
    Matrix C(nlin(),B.nlin());
#ifdef HAVE_BLAS
    DGEMM(CblasNoTrans,CblasTrans,
        (int)C.nlin(),(int)C.ncol(),(int)p,
        1.,data(),(int)nlin(),
        B.data(),(int)B.nlin(),
        0.,C.data(),(int)C.nlin());
#else
    for (size_t i=0;i<C.nlin();i++)
        for (size_t j=0;j<C.ncol();j++) {
            C(i,j)=0;
            for (size_t k=0;k<p;k++)
                C(i,j)+=(*this)(i,k)*B(j,k);
        }
#endif
        return C;
}

Matrix Matrix::tmultt(const Matrix &B) const
{
    assert(nlin()==B.ncol());
    size_t p=nlin();
    Matrix C(ncol(),B.nlin());
#ifdef HAVE_BLAS
    DGEMM(CblasTrans,CblasTrans,
        (int)C.nlin(),(int)C.ncol(),(int)p,
        1.,data(),(int)nlin(),
        B.data(),(int)B.nlin(),
        0.,C.data(),(int)C.nlin());
#else
    for (size_t i=0;i<C.nlin();i++)
        for (size_t j=0;j<C.ncol();j++) {
            C(i,j)=0;
            for (size_t k=0;k<p;k++)
                C(i,j)+=(*this)(k,i)*B(j,k);
        }
#endif
        return C;
}

Matrix Matrix::operator*(const SymMatrix &B) const
{
    assert(ncol()==B.ncol());
    Matrix C(nlin(),B.ncol());

#ifdef HAVE_BLAS
    Matrix D(B);
    DSYMM(CblasRight,  CblasUpper
        , (int)nlin(), (int)D.ncol(),
        1. , D.data(), (int)D.ncol(),
        data(), (int)nlin(),
        0, C.data(),(int)C.nlin());
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

Matrix Matrix::operator*(double x) const {
    Matrix C(nlin(),ncol());
    for (size_t k=0; k<nlin()*ncol(); k++) C.data()[k] = data()[k]*x;
    return C;
}

Matrix Matrix::operator/(double x) const {
    Matrix C(nlin(),ncol());
    for (size_t k=0; k<nlin()*ncol(); k++) C.data()[k] = data()[k]/x;
    return C;
}

Matrix Matrix::operator+(const Matrix &B) const
{
    assert(ncol()==B.ncol());
    assert(nlin()==B.nlin());
    Matrix C(*this,DEEP_COPY);
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(nlin()*ncol()), 1.0, B.data(), 1, C.data() , 1);
#else
    for (size_t i=0;i<nlin()*ncol();i++)
        C.data()[i]+=B.data()[i];
#endif
    return C;
}

Matrix Matrix::operator-(const Matrix &B) const
{
    assert(ncol()==B.ncol());
    assert(nlin()==B.nlin());
    Matrix C(*this,DEEP_COPY);
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(nlin()*ncol()), -1.0, B.data(), 1, C.data() , 1);
#else
    for (size_t i=0;i<nlin()*ncol();i++)
        C.data()[i]-=B.data()[i];
#endif
    return C;
}

void Matrix::operator+=(const Matrix &B)
{
    assert(ncol()==B.ncol());
    assert(nlin()==B.nlin());
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(nlin()*ncol()), 1.0, B.data(), 1, data() , 1);
#else
    for (size_t i=0;i<nlin()*ncol();i++)
        data()[i]+=B.data()[i];
#endif
}

void Matrix::operator-=(const Matrix &B)
{
    assert(ncol()==B.ncol());
    assert(nlin()==B.nlin());
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(nlin()*ncol()), -1.0, B.data(), 1, data() , 1);
#else
    for (size_t i=0;i<nlin()*ncol();i++)
        data()[i]-=B.data()[i];
#endif
}

void Matrix::operator*=(double x) {
    for (size_t k=0; k<nlin()*ncol(); k++) data()[k] *= x;
}

void Matrix::operator/=(double x) {
    for (size_t k=0; k<nlin()*ncol(); k++) data()[k] /= x;
}

double Matrix::dot(const Matrix& b) const {
    assert(nlin()==b.nlin()&&ncol()==b.ncol());
#ifdef HAVE_BLAS
    return BLAS(ddot,DDOT)((int)(nlin()*ncol()),data(),1,b.data(),1);
#else
    double s=0;
    for (size_t i=0;i<nlin()*ncol();i++)
        s+=data()[i]*b.data()[i];
    return s;
#endif
}

double Matrix::frobenius_norm() const {
#ifdef HAVE_LAPACK
    double info;
    Matrix b(*this,DEEP_COPY);
    return DLANGE('F',nlin(),ncol(),b.data(),nlin(),&info);
#else
    double d=0;
    for (size_t i=0; i<nlin()*ncol(); i++) d+=data()[i]*data()[i];
    return sqrt(d);
#endif
}

Vector Matrix::mean() const {
    Vector v(ncol()); v.set(0);
    for(size_t j = 0; j < ncol(); ++j) {
        for(size_t i = 0; i < nlin(); ++i) {
            v(j) += this->operator()(i,j);
        }
    }
    for(size_t j = 0; j < ncol(); ++j) {
        v(j) = v(j) / nlin();
    }
    return v;
}

Vector Matrix::tmean() const {
    Vector v(nlin()); v.set(0);
    for(size_t j = 0; j < ncol(); ++j) {
        for(size_t i = 0; i < nlin(); ++i) {
            v(i) += this->operator()(i,j);
        }
    }
    for(size_t i = 0; i < nlin(); ++i) {
        v(i) = v(i) / ncol();
    }
    return v;
}

void Matrix::info() const {
    if ((nlin() == 0) && (ncol() == 0)) {
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
        for(size_t j = 0; j < ncol(); ++j)
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
        for(size_t j = 0; j < std::min(ncol(),(size_t) 5); ++j)
        {
            std::cout << this->operator()(i,j) << " " ;
        }
        std::cout << std::endl ;
    }
}

// =======
// = IOs =
// =======

void Matrix::loadBin( const char *filename )
{
    maths::ifstream ifs(filename);
    ifs >> maths::format("binary") >> *this;
}

void Matrix::saveBin( const char *filename ) const
{
    maths::ofstream ofs(filename);
    ofs << maths::format("binary") << *this;
}

void Matrix::loadTxt( const char *filename )
{
    maths::ifstream ifs(filename);
    ifs >> maths::format("ascii") >> *this;
}

void Matrix::saveTxt( const char *filename ) const
{
    maths::ofstream ofs(filename);
    ofs << maths::format("ascii") << *this;
}

void Matrix::loadMat(const char *filename)
{
    maths::ifstream ifs(filename);
    ifs >> maths::format("matlab") >> *this;
}

void Matrix::saveMat( const char *filename ) const
{
    maths::ofstream ofs(filename);
    ofs << maths::format("matlab") << *this;
}

void Matrix::loadBrainvisa(const char *filename)
{
    maths::ifstream ifs(filename);
    ifs >> maths::format("tex") >> *this;
}

void Matrix::saveBrainvisa( const char *filename ) const
{
    maths::ofstream ofs(filename);
    ofs << maths::format("tex") << *this;
}

void Matrix::load( const char *filename ) {
    try {
        maths::ifstream ifs(filename);
        ifs >> *this;
    }
    catch (std::string s) {
        std::cout << s << std::endl;
    }
}

void Matrix::save( const char *filename ) const {
    try {
        maths::ofstream ofs(filename);
        ofs << *this;
    }
    catch (std::string s) {
        std::cout << s << std::endl;
    }
}
