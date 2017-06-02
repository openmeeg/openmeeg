/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

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

#pragma once

#include <OpenMEEGMathsConfig.h>
#include <iostream>
#include <cstdlib>
#include <string>

#include <linop.h>
#include <MathsIO.H>
#include <symmatrix.h>

namespace OpenMEEG {

    class SparseMatrix;
    class SymMatrix;
    class Vector;

    /** \brief  Matrix class

        Matrix class
    **/
    class OPENMEEGMATHS_EXPORT Matrix: public LinOp {
    protected:

        friend class Vector;

        utils::RCPtr<LinOpValue> value;

        explicit Matrix(const Matrix& A,const size_t M): LinOp(A.nlin(),M,FULL,2),value(A.value) { }

    public:

        Matrix(): LinOp(0,0,FULL,2),value() { }
        Matrix(const char* fname): LinOp(0,0,FULL,2),value() { this->load(fname); }
        Matrix(const size_t M,const size_t N): LinOp(M,N,FULL,2),value(new LinOpValue(N*M)) { }
        Matrix(const Matrix& A,const DeepCopy): LinOp(A.nlin(),A.ncol(),FULL,2),value(new LinOpValue(A.size(),A.data())) { }

        explicit Matrix(const SymMatrix& A);
        explicit Matrix(const SparseMatrix& A);

        Matrix(const Vector& v,const size_t M,const size_t N);

        void alloc_data()                       { value = new LinOpValue(size());      }
        void reference_data(const double* vals) { value = new LinOpValue(size(),vals); }

        /** \brief Test if Matrix is empty
            \return true if Matrix is empty
            \sa
        **/
        bool empty() const { return value->empty(); }

        /** \brief Get Matrix size
            \return number of values (nb lines x nb columns)
            \sa
        **/
        size_t size() const { return nlin()*ncol(); };

        /** \brief Get Matrix data
            \return pointer to Matrix values
            \sa
        **/
        double* data() const { return value->data; }

        /** \brief Get Matrix value
            \return value in Matrix
            \sa
        **/
        inline double operator()(size_t i,size_t j) const ;

        /** \brief Get Matrix value
            \return reference to value in Matrix
            \sa
        **/
        inline double& operator()(size_t i,size_t j) ;

        Matrix submat(size_t istart, size_t isize, size_t jstart, size_t jsize) const;
        void insertmat(size_t istart, size_t jstart, const Matrix& B);
        Vector getcol(size_t j) const;
        void setcol(size_t j, const Vector& v);
        Vector getlin(size_t i) const;
        void setlin(size_t i, const Vector& v);

        const Matrix& set(const double d);

        Matrix operator*(const Matrix& B) const;
        Matrix operator*(const SymMatrix& B) const;
        Matrix operator*(const SparseMatrix& B) const;
        Matrix operator+(const Matrix& B) const;
        Matrix operator-(const Matrix& B) const;
        Matrix operator*(double x) const;
        Matrix operator/(double x) const;
        void operator+=(const Matrix& B);
        void operator-=(const Matrix& B);
        void operator*=(double x);
        void operator/=(double x);

        Vector operator*(const Vector& v) const;
        Vector tmult(const Vector &v) const;
        Matrix tmult(const Matrix &m) const;
        Matrix multt(const Matrix &m) const;
        Matrix tmultt(const Matrix &m) const;

        Vector mean() const;
        Vector tmean() const;

        Matrix transpose() const;
        Matrix inverse() const;
        Matrix pinverse(double reltol=0) const;
        void svd(Matrix &U, SparseMatrix &S, Matrix &V, bool complete=true) const;

        /** \brief Get Matrix Frobenius norm
            \return norm value
            \sa
        **/
        double frobenius_norm() const;
        double dot(const Matrix& B) const;

        /** \brief Save Matrix to file (Format set using file name extension)
            \sa
        **/
        void save(const char *filename) const;

        /** \brief Load Matrix from file (Format set using file name extension)
            \sa
        **/
        void load(const char *filename);

        void save(const std::string& s) const { save(s.c_str()); }
        void load(const std::string& s)       { load(s.c_str()); }

        /** \brief Print info on Matrix
            \sa
        **/
        void info() const;

        friend class SparseMatrix;
        friend class SymMatrix;
    };

    inline double Matrix::operator()(size_t i,size_t j) const {
        om_assert(i<nlin() && j<ncol());
        return value->data[i+nlin()*j];
    }
    inline double& Matrix::operator()(size_t i,size_t j) {
        om_assert(i<nlin() && j<ncol());
        return value->data[i+nlin()*j];
    }
    
    inline double Matrix::frobenius_norm() const {
    #ifdef HAVE_LAPACK
    if ( nlin()*ncol() != 0 ) {
        double work;
        return DLANGE('F', sizet_to_int(nlin()), sizet_to_int(ncol()), data(), sizet_to_int(nlin()), &work);
    } else {
        return 0;
    }
    #else
        double d = 0.;
        for (size_t i=0; i<nlin()*ncol(); i++) d+=data()[i]*data()[i];
        return sqrt(d);
    #endif
    }

    inline Vector Matrix::operator*(const Vector &v) const {
        om_assert(ncol()==v.nlin());
        Vector y(nlin());
    #ifdef HAVE_BLAS
        DGEMV(CblasNoTrans,sizet_to_int(nlin()),sizet_to_int(ncol()),1.0,data(),sizet_to_int(nlin()),v.data(),1,0.,y.data(),1);
    #else
        for (size_t i=0;i<nlin();i++) {
            y(i) = 0;
            for (size_t j=0;j<ncol();j++)
                y(i) += (*this)(i,j)*v(j);
        }
    #endif

        return y;
    }

    inline Matrix Matrix::submat(size_t istart, size_t isize, size_t jstart, size_t jsize) const {
        om_assert (istart+isize<=nlin() && jstart+jsize<=ncol());

        Matrix a(isize,jsize);
        const int sz = sizet_to_int(isize);
        
        for (size_t j=0; j<jsize; j++) {
    #ifdef HAVE_BLAS
            BLAS(dcopy,DCOPY)(sz,data()+istart+(jstart+j)*nlin(),1,a.data()+j*isize,1);
    #else
            for (size_t i=0; i<isize; i++)
                a(i,j) = (*this)(istart+i,jstart+j);
    #endif
        }
        return a;
    }

    inline void Matrix::insertmat(size_t istart, size_t jstart, const Matrix& B) {
        om_assert (istart+B.nlin()<=nlin() && jstart+B.ncol()<=ncol() );
        for (size_t j=0; j<B.ncol(); j++) {
            for (size_t i=0; i<B.nlin(); i++) {
                (*this)(istart+i,jstart+j)=B(i,j);
            }
        }
    }

    inline Vector Matrix::getcol(size_t j) const {
        om_assert(j<ncol());
        Vector v(nlin());
    #ifdef HAVE_BLAS
        BLAS(dcopy,DCOPY)(sizet_to_int(nlin()),data()+nlin()*j,1,v.data(),1);
    #else
        for (size_t i=0;i<nlin();i++) v.data()[i]=data()[i+nlin()*j];
    #endif
        return v;
    }

    inline Vector Matrix::getlin(size_t i) const {
        om_assert(i<nlin());
        Vector v(ncol());
    #ifdef HAVE_BLAS
        BLAS(dcopy,DCOPY)(sizet_to_int(ncol()),data()+i,sizet_to_int(nlin()),v.data(),1);
    #else
        for (size_t j=0;j<ncol();j++) v.data()[j]=data()[i+nlin()*j];
    #endif
        return v;
    }

    inline void Matrix::setcol(size_t j,const Vector& v) {
        om_assert(v.size()==nlin() && j<ncol());
    #ifdef HAVE_BLAS
        BLAS(dcopy,DCOPY)(sizet_to_int(nlin()),v.data(),1,data()+nlin()*j,1);
    #else
        for (size_t i=0;i<nlin();i++) data()[i+nlin()*j]=v.data()[i];
    #endif
    }

    inline void Matrix::setlin(size_t i,const Vector& v) {
        om_assert(v.size()==ncol() && i<nlin());
    #ifdef HAVE_BLAS
        BLAS(dcopy,DCOPY)(sizet_to_int(ncol()),v.data(),1,data()+i,sizet_to_int(nlin()));
    #else
        for (size_t j=0;j<ncol();j++) data()[i+nlin()*j]=v.data()[j];
    #endif
    }

    inline Vector Matrix::tmult(const Vector &v) const {
        om_assert(nlin()==v.nlin());
        Vector y(ncol());
    #ifdef HAVE_BLAS
        DGEMV(CblasTrans,sizet_to_int(nlin()),sizet_to_int(ncol()),1.,data(),sizet_to_int(nlin()),v.data(),1,0.,y.data(),1);
    #else
        for (size_t i=0;i<ncol();i++) {
            y(i)=0;
            for (size_t j=0;j<nlin();j++)
                y(i)+=(*this)(j,i)*v(j);
        }
    #endif

        return y;
    }

    inline Matrix Matrix::inverse() const {
    #ifdef HAVE_LAPACK
        om_assert(nlin()==ncol());
        Matrix invA(*this,DEEP_COPY);
        // LU
        #if defined(CLAPACK_INTERFACE)
            #if defined(__APPLE__) && defined(USE_VECLIB) // Apple Veclib Framework (Handles 32 and 64 Bits)
                __CLPK_integer *pivots = new __CLPK_integer[ncol()];
                __CLPK_integer Info = 0;
                __CLPK_integer nlin_local = invA.nlin();
                __CLPK_integer nlin_local2 = invA.nlin();
                __CLPK_integer ncol_local = invA.ncol();
                __CLPK_integer sz = invA.ncol()*64;
                DGETRF(nlin_local,ncol_local,invA.data(),nlin_local2,pivots,Info);
                double *work=new double[sz];
                DGETRI(ncol_local,invA.data(),ncol_local,pivots,work,sz,Info);
                delete[] pivots;
                delete[] work;
                om_assert(Info==0);
            #else
                BLAS_INT *pivots=new BLAS_INT[sizet_to_int(ncol())];
                DGETRF(sizet_to_int(invA.nlin()),sizet_to_int(invA.ncol()),invA.data(),sizet_to_int(invA.nlin()),pivots);
                DGETRI(sizet_to_int(invA.ncol()),invA.data(),sizet_to_int(invA.ncol()),pivots);
                delete[] pivots;
            #endif
        #else
            int Info = 0;
            int *pivots=new int[sizet_to_int(ncol())];
            DGETRF(sizet_to_int(invA.nlin()),sizet_to_int(invA.ncol()),invA.data(),sizet_to_int(invA.nlin()),pivots,Info);
            const unsigned sz = invA.ncol()*64;
            double *work=new double[sz];
            DGETRI(sizet_to_int(invA.ncol()),invA.data(),sizet_to_int(invA.ncol()),pivots,work,sz,Info);
            delete[] pivots;
            delete[] work;
            om_assert(Info==0);
        #endif
        return invA;
    #else
        std::cerr << "!!!!! Inverse not implemented !!!!!" << std::endl;
        exit(1);
    #endif
    }

    inline Matrix Matrix::operator *(const Matrix &B) const {
        om_assert(ncol()==B.nlin());
        size_t p=ncol();
        Matrix C(nlin(),B.ncol());
    #ifdef HAVE_BLAS
        DGEMM(CblasNoTrans,CblasNoTrans,
            sizet_to_int(C.nlin()),sizet_to_int(C.ncol()),sizet_to_int(p),
            1.,data(),sizet_to_int(nlin()),
            B.data(),sizet_to_int(B.nlin()),
            0.,C.data(),sizet_to_int(C.nlin()));
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
    
    inline Matrix Matrix::tmult(const Matrix &B) const {
        om_assert(nlin()==B.nlin());
        size_t p=nlin();
        Matrix C(ncol(),B.ncol());
    #ifdef HAVE_BLAS
        DGEMM(CblasTrans,CblasNoTrans,
            sizet_to_int(C.nlin()),sizet_to_int(C.ncol()),sizet_to_int(p),
            1.,data(),sizet_to_int(nlin()),
            B.data(),sizet_to_int(B.nlin()),
            0.,C.data(),sizet_to_int(C.nlin()));
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

    inline Matrix Matrix::multt(const Matrix &B) const {
        om_assert(ncol()==B.ncol());
        size_t p=ncol();
        Matrix C(nlin(),B.nlin());
    #ifdef HAVE_BLAS
        DGEMM(CblasNoTrans,CblasTrans,
            sizet_to_int(C.nlin()),sizet_to_int(C.ncol()),sizet_to_int(p),
            1.,data(),sizet_to_int(nlin()),
            B.data(),sizet_to_int(B.nlin()),
            0.,C.data(),sizet_to_int(C.nlin()));
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

    inline Matrix Matrix::tmultt(const Matrix &B) const {
        om_assert(nlin()==B.ncol());
        size_t p=nlin();
        Matrix C(ncol(),B.nlin());
    #ifdef HAVE_BLAS
        DGEMM(CblasTrans,CblasTrans,
            sizet_to_int(C.nlin()),sizet_to_int(C.ncol()),sizet_to_int(p),
            1.,data(),sizet_to_int(nlin()),
            B.data(),sizet_to_int(B.nlin()),
            0.,C.data(),sizet_to_int(C.nlin()));
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

    inline Matrix Matrix::operator*(const SymMatrix &B) const {
        om_assert(ncol()==B.ncol());
        Matrix C(nlin(),B.ncol());

    #ifdef HAVE_BLAS
        Matrix D(B);
        const int n = sizet_to_int(nlin());
        const int m = sizet_to_int(D.ncol());
        const int l = sizet_to_int(C.nlin());
        DSYMM(CblasRight,CblasUpper,n,m,1.,D.data(),m,data(),n,0.,C.data(),l);
    #else
        for (size_t j=0;j<B.ncol();j++)
            for (size_t i=0;i<ncol();i++) {
                C(i,j)=0;
                for (size_t k=0;k<ncol();k++)
                    C(i,j)+=(*this)(i,k)*B(k,j);
            }
    #endif
            return C;
    }
    
    inline Matrix Matrix::operator+(const Matrix &B) const {
        om_assert(ncol()==B.ncol());
        om_assert(nlin()==B.nlin());
        Matrix C(*this,DEEP_COPY);
    #ifdef HAVE_BLAS
        BLAS(daxpy,DAXPY)(sizet_to_int(nlin()*ncol()), 1.0, B.data(), 1, C.data() , 1);
    #else
        for (size_t i=0;i<nlin()*ncol();i++)
            C.data()[i]+=B.data()[i];
    #endif
        return C;
    }

    inline Matrix Matrix::operator-(const Matrix &B) const {
        om_assert(ncol()==B.ncol());
        om_assert(nlin()==B.nlin());
        Matrix C(*this,DEEP_COPY);
    #ifdef HAVE_BLAS
        BLAS(daxpy,DAXPY)(sizet_to_int(nlin()*ncol()), -1.0, B.data(), 1, C.data() , 1);
    #else
        for (size_t i=0;i<nlin()*ncol();i++)
            C.data()[i]-=B.data()[i];
    #endif
        return C;
    }

    inline void Matrix::operator+=(const Matrix &B) {
        om_assert(ncol()==B.ncol());
        om_assert(nlin()==B.nlin());
    #ifdef HAVE_BLAS
        BLAS(daxpy,DAXPY)(sizet_to_int(nlin()*ncol()), 1.0, B.data(), 1, data() , 1);
    #else
        for (size_t i=0;i<nlin()*ncol();i++)
            data()[i]+=B.data()[i];
    #endif
    }

    inline void Matrix::operator-=(const Matrix &B) {
        om_assert(ncol()==B.ncol());
        om_assert(nlin()==B.nlin());
    #ifdef HAVE_BLAS
        BLAS(daxpy,DAXPY)(sizet_to_int(nlin()*ncol()), -1.0, B.data(), 1, data() , 1);
    #else
        for (size_t i=0;i<nlin()*ncol();i++)
            data()[i]-=B.data()[i];
    #endif
    }
    
    inline double Matrix::dot(const Matrix& b) const {
        om_assert(nlin()==b.nlin()&&ncol()==b.ncol());
    #ifdef HAVE_BLAS
        return BLAS(ddot,DDOT)(sizet_to_int(nlin()*ncol()),data(),1,b.data(),1);
    #else
        double s=0;
        for (size_t i=0;i<nlin()*ncol();i++)
            s+=data()[i]*b.data()[i];
        return s;
    #endif
    }
}
