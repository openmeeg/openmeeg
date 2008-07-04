/* FILE: $Id$ */

/*
Project Name : OpenMEEG

author            : $Author$
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
#include "matrice_dcl.h"
#include "symmatrice_dcl.h"
#include "sparse_matrice_dcl.h"
#include "vecteur_dcl.h"
#include "sparse_matrice_dcl.h"

matrice::matrice() : MatrixBase(0,0,FULL,TWO),t(0),count(0) { }
matrice::matrice(const char* fname) : MatrixBase(0,0,FULL,TWO),t(0),count(0) {
    this->load(fname);
}

matrice::matrice(size_t M,size_t N) : MatrixBase(M,N,FULL,TWO),t(0),count(0) { alloc_data(); }
matrice::matrice(double* T, int* COUNT, size_t M, size_t N) : MatrixBase(M,N,FULL,TWO),t(T),count(COUNT) {(*count)++;}

bool matrice::empty() const { return t==0;}
double* matrice::data() const {return t;}
int* matrice::DangerousGetCount () const {return count;}

matrice matrice::duplicate() const
{
    matrice A;
    if (t) {
        A.nlin() = nlin();
        A.ncol() = ncol();
        A.alloc_data();
        copyout(A.t);
    }
    return A;
}

void matrice::copyin(const matrice& A)
{
    if (t) {
        assert(nlin()==A.nlin() && ncol()==A.ncol());
        copyin(A.t);
    }
}

// double matrice::operator[](size_t i) const // TODO : remove?
// {
//     assert(i<nlin()*ncol());
//     return t[i];
// }
// double& matrice::operator[](size_t i)
// {
//     assert(i<nlin()*ncol());
//     return t[i];
// }

void matrice::alloc_data() {
    if (t!=0)
        destroy();
    t = new double[size()];
    count = new int[1];
    (*count) = 1;
}

void matrice::destroy()
{
    if (t!=0) {
        (*count)--;
        if ((*count)==0) {
            delete[] t;
            delete[] count;
        }
    }
}

void matrice::copy(const matrice& A)
{
    t=A.t;
    nlin()=A.nlin();
    ncol()=A.ncol();
    if (t) {
        count = A.count;
        (*count)++;
    }
}

void matrice::copyout(double * p) const {
    if (!t) return;
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)(nlin()*ncol()),t,1,p,1);
#else
    for (size_t i=0;i<nlin()*ncol();i++)
        p[i]=t[i];
#endif
}

void matrice::copyin(const double * p) {
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)(nlin()*ncol()),p,1,t,1);
#else
    for (size_t i=0;i<nlin()*ncol();i++)
        t[i]=p[i];
#endif
}

const matrice& matrice::operator=(const matrice& A)
{
    destroy();
    copy(A);
    return *this;
}

const matrice& matrice::set(const double d)
{
    for(size_t i=0;i<ncol()*nlin();i++) t[i]=d;
    return *this;
}

matrice::matrice(const matrice& A): MatrixBase(A.nlin(),A.ncol(),FULL,TWO)
{
    copy(A);
}

matrice::matrice(const symmatrice& A): MatrixBase(A.nlin(),A.ncol(),FULL,TWO)
{
    nlin() = A.nlin();
    ncol() = A.ncol();
    alloc_data();
    for (size_t j=0; j<ncol(); j++)
        for (size_t i=0; i<nlin(); i++)
            (*this)(i,j)=A(i,j);
}

matrice::matrice(const vecteur& v, size_t M, size_t N): MatrixBase(M,N,FULL,TWO)
{
    assert(M*N==v.size());
    t=v.data();
    if (t) {
        count=v.DangerousGetCount();
        (*count)++;
    }
}

vecteur matrice::operator *(const vecteur &v) const
{
    assert(ncol()==v.nlin());
    vecteur y(nlin());
#ifdef HAVE_BLAS
    DGEMV(CblasNoTrans,(int)nlin(),(int)ncol(),1.0,t,(int)nlin(),v.t,1,0.,y.t,1);
#else
    for (size_t i=0;i<nlin();i++) {
        y(i)=0;
        for (size_t j=0;j<ncol();j++)
            y(i)+=(*this)(i,j)*v(j);
    }
#endif

    return y;
}

matrice matrice::submat(size_t istart, size_t isize, size_t jstart, size_t jsize) const {
    assert (istart+isize<=nlin() && jstart+jsize<=ncol());

    matrice a(isize,jsize);
    for (size_t j=0; j<jsize; j++)
#ifdef HAVE_BLAS
        BLAS(dcopy,DCOPY)((int)(isize),t+istart+(jstart+j)*nlin(),1,a.t+j*isize,1);
#elif USE_ACML
        dcopy((int)(isize),t+istart+(jstart+j)*nlin(),1,a.t+j*isize,1);
#else
        for (size_t i=0; i<isize; i++)
            a(i,j)=(*this)(istart+i,jstart+j);
#endif
    return a;
}

vecteur matrice::getcol(size_t j) const {
    assert(j<ncol());
    vecteur v(nlin());
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)nlin(),t+nlin()*j,1,v.t,1);
#else
    for (size_t i=0;i<nlin();i++) v.t[i]=t[i+nlin()*j];
#endif
    return v;
}

vecteur matrice::getlin(size_t i) const {
    assert(i<nlin());
    vecteur v(ncol());
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)ncol(),t+i,(int)nlin(),v.t,1);
#else
    for (size_t j=0;j<ncol();j++) v.t[j]=t[i+nlin()*j];
#endif
    return v;
}

void matrice::setcol(size_t j, const vecteur& v) {
    assert(v.size()==nlin() && j<ncol());
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)nlin(),v.t,1,t+nlin()*j,1);
#else
    for (size_t i=0;i<nlin();i++) t[i+nlin()*j]=v.t[i];
#endif
}

void matrice::setlin(size_t i, const vecteur& v) {
    assert(v.size()==ncol() && i<nlin());
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)ncol(),v.t,1,t+i,(int)nlin());
#else
    for (size_t j=0;j<ncol();j++) t[i+nlin()*j]=v.t[j];
#endif
}

vecteur matrice::tmult(const vecteur &v) const
{
    assert(nlin()==v.nlin());
    vecteur y(ncol());
#ifdef HAVE_BLAS
    DGEMV(CblasTrans,(int)nlin(),(int)ncol(),1.,t,(int)nlin(),v.t,1,0.,y.t,1);
#else
    for (size_t i=0;i<ncol();i++) {
        y(i)=0;
        for (size_t j=0;j<nlin();j++)
            y(i)+=(*this)(j,i)*v(j);
    }
#endif

    return y;
}

matrice matrice::inverse() const
{
#ifdef HAVE_LAPACK
    assert(nlin()==ncol());
    matrice invA=duplicate();
    // LU
    int *pivots=new int[ncol()];
    int info;
    DGETRF(invA.nlin(),invA.ncol(),invA.t,invA.nlin(),pivots,info);
    // Inverse
    int size=(int)invA.ncol()*64;
    double *work=new double[size];
    DGETRI(invA.ncol(),invA.t,invA.ncol(),pivots,work,size,info);
    delete[] pivots;
    delete[] work;
    return invA;
#else
    std::cerr << "!!!!! Inverse not implemented !!!!!" << std::endl;
    exit(1);
#endif
}

matrice matrice::pinverse(double tolrel) const {
#if defined(HAVE_BLAS) && defined(HAVE_LAPACK)
    if(ncol() > nlin()) return transpose().pinverse().transpose();
    else {
        matrice result(ncol(),nlin());
        matrice U,S,V;
        svd(U,S,V);
        double maxs=0;
        int mimi=(int)std::min(S.nlin(),S.ncol());
        for(int i=0;i<mimi;i++) maxs=std::max(S(i,i),maxs);
        if (tolrel==0) tolrel=DBL_EPSILON;
        double tol = std::max(nlin(),ncol()) * maxs * tolrel;
        int r=0; for(int i=0;i<mimi;i++) if(S(i,i)>tol) r++;
        if (r == 0)
        {
            result.set(0.);
            return result;
        }
        else
        {
            matrice s(r,r); s.set(0);
            for(int i=0;i<r;i++) s(i,i)=1.0/S(i,i);
            matrice Vbis;
            Vbis.DangerousBuild(V.t,V.nlin(),r);
            matrice Ubis;
            Ubis.DangerousBuild(U.t,U.nlin(),r);
            result=Vbis*s*Ubis.transpose();
            Vbis.DangerousKill();
            Ubis.DangerousKill();
            return result;
        }
    }
#else
    std::cerr << "pinv not implemented without blas/lapack" << std::endl;
    exit(1);
#endif
}

matrice matrice::transpose() const {
    matrice result(ncol(),nlin());
    for(size_t i=0;i<nlin();i++) for(size_t j=0;j<ncol();j++) result(j,i)=(*this)(i,j);
    return result;
}

void matrice::svd(matrice &U,matrice &S, matrice &V) const {
#ifdef HAVE_LAPACK
    matrice cpy=duplicate();
    int mimi=(int)std::min(nlin(),ncol());
    U=matrice(nlin(),ncol()); U.set(0);
    V=matrice(ncol(),ncol()); V.set(0);
    S=matrice(ncol(),ncol()); S.set(0);
    double *s=new double[mimi];
    int lwork=4 *mimi*mimi + (int)std::max(nlin(),ncol()) + 9*mimi;
    double *work=new double[lwork];
    int *iwork=new int[8*mimi];
    int info;
    DGESDD('S',nlin(),ncol(),cpy.t,nlin(),s,U.t,U.nlin(),V.t,V.nlin(),work,lwork,iwork,info);
    for(int i=0;i<mimi;i++) S(i,i)=s[i];
    V=V.transpose();
    delete[] s;
    delete[] work;
    delete[] iwork;
#else
    std::cerr<<"svd not implemented without blas/lapack"<<std::endl;
#endif
}

matrice matrice::operator *(const matrice &B) const
{
    assert(ncol()==B.nlin());
    size_t p=ncol();
    matrice C(nlin(),B.ncol());
#ifdef HAVE_BLAS
    DGEMM(CblasNoTrans,CblasNoTrans,
        (int)C.nlin(),(int)C.ncol(),(int)p,
        1.,t,(int)nlin(),
        B.t,(int)B.nlin(),
        0.,C.t,(int)C.nlin());
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

matrice matrice::operator *(const sparse_matrice &mat) const
{
    assert(ncol()==mat.nlin());
    matrice out(nlin(),mat.ncol());
    out.set(0.0);

    sparse_matrice::const_iterator it;
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

matrice matrice::tmult(const matrice &B) const
{
    assert(nlin()==B.nlin());
    size_t p=nlin();
    matrice C(ncol(),B.ncol());
#ifdef HAVE_BLAS
    DGEMM(CblasTrans,CblasNoTrans,
        (int)C.nlin(),(int)C.ncol(),(int)p,
        1.,t,(int)nlin(),
        B.t,(int)B.nlin(),
        0.,C.t,(int)C.nlin());
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

matrice matrice::multt(const matrice &B) const
{
    assert(ncol()==B.ncol());
    size_t p=ncol();
    matrice C(nlin(),B.nlin());
#ifdef HAVE_BLAS
    DGEMM(CblasNoTrans,CblasTrans,
        (int)C.nlin(),(int)C.ncol(),(int)p,
        1.,t,(int)nlin(),
        B.t,(int)B.nlin(),
        0.,C.t,(int)C.nlin());
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

matrice matrice::tmultt(const matrice &B) const
{
    assert(nlin()==B.ncol());
    size_t p=nlin();
    matrice C(ncol(),B.nlin());
#ifdef HAVE_BLAS
    DGEMM(CblasTrans,CblasTrans,
        (int)C.nlin(),(int)C.ncol(),(int)p,
        1.,t,(int)nlin(),
        B.t,(int)B.nlin(),
        0.,C.t,(int)C.nlin());
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

matrice matrice::operator *(const symmatrice &B) const
{
    assert(ncol()==B.ncol());
    matrice C(nlin(),B.ncol());

#ifdef HAVE_BLAS
    matrice D(B);
    DSYMM(CblasRight,  CblasUpper
        , (int)nlin(), (int)D.ncol(),
        1. , D.t, (int)D.ncol(),
        t, (int)nlin(),
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

matrice matrice::operator *(double x) const {
    matrice C(nlin(),ncol());
    for (size_t k=0; k<nlin()*ncol(); k++) C.t[k] = t[k]*x;
    return C;
}

matrice matrice::operator /(double x) const {
    matrice C(nlin(),ncol());
    for (size_t k=0; k<nlin()*ncol(); k++) C.t[k] = t[k]/x;
    return C;
}

matrice matrice::operator +(const matrice &B) const
{
    assert(ncol()==B.ncol());
    assert(nlin()==B.nlin());
    matrice C=duplicate();
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(nlin()*ncol()), 1.0, B.t, 1, C.t , 1);
#else
    for (size_t i=0;i<nlin()*ncol();i++)
        C.t[i]+=B.t[i];
#endif
    return C;
}

matrice matrice::operator -(const matrice &B) const
{
    assert(ncol()==B.ncol());
    assert(nlin()==B.nlin());
    matrice C=duplicate();
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(nlin()*ncol()), -1.0, B.t, 1, C.t , 1);
#else
    for (size_t i=0;i<nlin()*ncol();i++)
        C.t[i]-=B.t[i];
#endif
    return C;
}

void matrice::operator +=(const matrice &B)
{
    assert(ncol()==B.ncol());
    assert(nlin()==B.nlin());
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(nlin()*ncol()), 1.0, B.t, 1, t , 1);
#else
    for (size_t i=0;i<nlin()*ncol();i++)
        t[i]+=B.t[i];
#endif
}

void matrice::operator -=(const matrice &B)
{
    assert(ncol()==B.ncol());
    assert(nlin()==B.nlin());
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(nlin()*ncol()), -1.0, B.t, 1, t , 1);
#else
    for (size_t i=0;i<nlin()*ncol();i++)
        t[i]-=B.t[i];
#endif
}

void matrice::operator *=(double x) {
    for (size_t k=0; k<nlin()*ncol(); k++) t[k] *= x;
}

void matrice::operator /=(double x) {
    for (size_t k=0; k<nlin()*ncol(); k++) t[k] /= x;
}

void matrice::DangerousBuild( double *pt, size_t i, size_t j)
{
    t=pt;
    nlin()=i;
    ncol()=j;
    count=new int[1];
    *count=1;
}

void matrice::DangerousKill ()
{
    delete[] count;
    t=0;
}

double matrice::dot(const matrice& b) const {
    assert(nlin()==b.nlin()&&ncol()==b.ncol());
#ifdef HAVE_BLAS
    return BLAS(ddot,DDOT)((int)(nlin()*ncol()),t,1,b.t,1);
#else
    double s=0;
    for (size_t i=0;i<nlin()*ncol();i++)
        s+=t[i]*b.t[i];
    return s;
#endif
}

double matrice::frobenius_norm() const {
#ifdef HAVE_LAPACK
    double info;
    matrice b=duplicate();
    return DLANGE('F',nlin(),ncol(),b.t,nlin(),&info);
#else
    double d=0;
    for (size_t i=0; i<nlin()*ncol(); i++) d+=t[i]*t[i];
    return sqrt(d);
#endif
}

vecteur matrice::mean() const {
    vecteur v(ncol()); v.set(0);
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

vecteur matrice::tmean() const {
    vecteur v(nlin()); v.set(0);
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

void matrice::info() const {
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

void matrice::loadBin( const char *filename )
{
    Maths::ifstream ifs(filename);
    ifs >> Maths::format("old_binary") >> *this;
}

void matrice::saveBin( const char *filename ) const
{
    Maths::ofstream ofs(filename);
    ofs << Maths::format("old_binary") << *this;
}

void matrice::loadTxt( const char *filename )
{
    Maths::ifstream ifs(filename);
    ifs >> Maths::format("ascii") >> *this;
}

void matrice::saveTxt( const char *filename ) const
{
    Maths::ofstream ofs(filename);
    ofs << Maths::format("ascii") << *this;
}

void matrice::loadMat(const char *filename)
{
    Maths::ifstream ifs(filename);
    ifs >> Maths::format("matlab") >> *this;
}

void matrice::saveMat( const char *filename ) const
{
    Maths::ofstream ofs(filename);
    ofs << Maths::format("matlab") << *this;
}

void matrice::load( const char *filename ) {
    try {
        Maths::ifstream ifs(filename);
        ifs >> *this;
    }
    catch (std::string s) {
        std::cout << s << std::endl;
    }
}

void matrice::save( const char *filename ) const {
    try {
        Maths::ofstream ofs(filename);
        ofs << *this;
    }
    catch (std::string s) {
        std::cout << s << std::endl;
    }
}

std::ostream& operator<<(std::ostream& f,const matrice &M) {
    for (size_t i=0;i<M.nlin();i++) {
        for (size_t j=0;j<M.ncol();j++) {
            f << M(i,j) << " ";
        }
        f << std::endl;
    }
    return f;
}
