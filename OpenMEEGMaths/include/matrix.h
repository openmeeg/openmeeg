// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

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

    /// \brief  Matrix class
    /// Matrix class

    class OPENMEEGMATHS_EXPORT Matrix: public LinOp {
    protected:

        friend class Vector;

        LinOpValue value;

        explicit Matrix(const Matrix& A,const Dimension M): LinOp(A.nlin(),M,FULL,2),value(A.value) { }

    public:

        Matrix(): LinOp(0,0,FULL,2),value() { }
        Matrix(const char* fname): LinOp(0,0,FULL,2),value() { load(fname); }
        Matrix(const std::string& fname): Matrix(fname.c_str()) { }
        Matrix(const Dimension M,const Dimension N): LinOp(M,N,FULL,2),value(N*M) { }
        Matrix(const Matrix& A,const DeepCopy): LinOp(A.nlin(),A.ncol(),FULL,2),value(A.size(),A.data()) { }

        explicit Matrix(const SymMatrix& A);
        explicit Matrix(const SparseMatrix& A);

        Matrix(const Vector& v,const Dimension M,const Dimension N);

        void alloc_data()                       { value = LinOpValue(size());      }
        void reference_data(const double* vals) { value = LinOpValue(size(),vals); }

        /// \brief Test if Matrix is empty
        /// \return true if Matrix is empty

        bool empty() const { return value.empty(); }

        /// \brief Get Matrix size
        /// \return number of values (nb lines x nb columns)

        size_t size() const { return static_cast<std::streamoff>(nlin())*static_cast<std::streamoff>(ncol()); };

        /// \brief Get Matrix data
        /// \return pointer to Matrix values

        double* data() const { return value.get(); }

    #if defined(SWIGPYTHON) || defined(WIN32)

        // This method is only needed for swig wrapping. But windows insists on having it
        // in the library despite the fact it is inlined.

        /// \brief Get a shared pointer on the Matrix data.
        
        std::shared_ptr<double[]> get_shared_data_ptr() { return value; }
    #endif

        /// \brief Get Matrix value
        /// \return value in Matrix

        double operator()(const Index i,const Index j) const {
            om_assert(i<nlin() && j<ncol());
            return value[i+nlin()*j];
        }

        /// \brief Get Matrix value
        /// \return reference to value in Matrix

        double& operator()(const Index i,const Index j) {
            om_assert(i<nlin() && j<ncol());
            return value[i+nlin()*j];
        }

        Matrix submat(const Index istart,const Index isize,const Index jstart,const Index jsize) const;
        void   insertmat(const Index istart,const Index jstart,const Matrix& B);

        Vector getcol(const Index j) const;
        void   setcol(const Index j,const Vector& v);

        Vector getlin(const Index i) const;
        void   setlin(const Index i,const Vector& v);

        const Matrix& set(const double d);

        Matrix operator*(const Matrix& B) const;
        Matrix operator*(const SymMatrix& B) const;
        Matrix operator*(const SparseMatrix& B) const;

        Matrix operator+(const Matrix& B) const {
            Matrix C(*this,DEEP_COPY);
            C += B;
            return C;
        }

        Matrix operator-(const Matrix& B) const {
            Matrix C(*this,DEEP_COPY);
            C -= B;
            return C;
        }

        Matrix operator*(double x) const;
        Matrix operator/(double x) const;
        inline void operator+=(const Matrix& B);
        inline void operator-=(const Matrix& B);
        void operator*=(double x);
        void operator/=(double x);

        Vector operator*(const Vector& v) const;
        Vector tmult(const Vector& v) const;
        Matrix tmult(const Matrix& m) const;
        Matrix multt(const Matrix& m) const;
        Matrix tmultt(const Matrix& m) const;

        Matrix transpose() const;
        Matrix inverse() const;
        Matrix pinverse(const double reltol=0.0) const;
        void svd(Matrix& U,SparseMatrix& S,Matrix& V,const bool complete=true) const;

        /// \brief Get Matrix Frobenius norm
        /// \return norm value

        double frobenius_norm() const;
        double dot(const Matrix& B) const;

        /// \brief Save Matrix to file (Format set using file name extension)

        void save(const char* filename) const;

        /// \brief Load Matrix from file (Format set using file name extension)

        void load(const char* filename);

        void save(const std::string& s) const { save(s.c_str()); }
        void load(const std::string& s)       { load(s.c_str()); }

        /// \brief Print info on Matrix

        void info() const;

        friend class SparseMatrix;
        friend class SymMatrix;
    };

    inline std::ostream& operator<<(std::ostream& os,const Matrix& M) {
        for (unsigned i=0; i<M.nlin(); ++i) {
            for (unsigned j=0; j<M.ncol(); ++j)
                os << M(i,j) << ' ';
            os << std::endl;
        }
        return os;
    }

    inline double Matrix::frobenius_norm() const {
        const size_t sz = size();
        if (sz==0)
            return 0.0;

    #ifdef HAVE_LAPACK
        double work;
        const BLAS_INT M = sizet_to_int(nlin());
        const BLAS_INT N = sizet_to_int(ncol());
        return DLANGE('F',M,N,data(),M,&work);
    #else
        double d = 0.0;
        for (size_t i=0; i<sz; i++)
            d += data()[i]*data()[i];
        return sqrt(d);
    #endif
    }

    inline Vector Matrix::operator*(const Vector& v) const {
        om_assert(ncol()==v.nlin());
        Vector res(nlin());
    #ifdef HAVE_BLAS
        const BLAS_INT M = sizet_to_int(nlin());
        const BLAS_INT N = sizet_to_int(ncol());
        DGEMV(CblasNoTrans,M,N,1.0,data(),M,v.data(),1,0.0,res.data(),1);
    #else
        res.set(0.0);
        for (Index j=0; j<ncol(); ++j)
            for (Index i=0; i<nlin(); ++i)
                res(i) += (*this)(i,j)*v(j);
    #endif

        return res;
    }

    inline Matrix Matrix::submat(const Index istart,const Index isize,const Index jstart,const Index jsize) const {
        om_assert(istart+isize<=nlin() && jstart+jsize<=ncol());

        Matrix res(isize,jsize);

        for (Index j=0; j<jsize; ++j) {
    #ifdef HAVE_BLAS
            const BLAS_INT sz = sizet_to_int(isize);
            BLAS(dcopy,DCOPY)(sz,data()+istart+(jstart+j)*nlin(),1,res.data()+j*isize,1);
    #else
            for (Index i=0; i<isize; ++i)
                res(i,j) = (*this)(istart+i,jstart+j);
    #endif
        }
        return res;
    }

    inline void Matrix::insertmat(const Index istart,const Index jstart,const Matrix& B) {
        om_assert(istart+B.nlin()<=nlin() && jstart+B.ncol()<=ncol() );

        for (Index j=0; j<B.ncol(); ++j)
            for (Index i=0; i<B.nlin(); ++i)
                (*this)(istart+i,jstart+j) = B(i,j);
    }

    inline Vector Matrix::getcol(const Index j) const {
        om_assert(j<ncol( ));
        Vector res(nlin());
    #ifdef HAVE_BLAS
        const BLAS_INT M = sizet_to_int(nlin());
        BLAS(dcopy,DCOPY)(M,data()+nlin()*j,1,res.data(),1);
    #else
        for (Index i=0; i<nlin(); ++i)
            res(i) = (*this)(i,j);
    #endif
        return res;
    }

    inline Vector Matrix::getlin(const Index i) const {
        om_assert(i<nlin());
        Vector res(ncol());
    #ifdef HAVE_BLAS
        const BLAS_INT M = sizet_to_int(nlin());
        const BLAS_INT N = sizet_to_int(ncol());
        BLAS(dcopy,DCOPY)(N,data()+i,M,res.data(),1);
    #else
        for (Index j=0; j<ncol(); ++j)
            res(j) = (*this)(i,j);
    #endif
        return res;
    }

    inline void Matrix::setcol(const Index j,const Vector& v) {
        om_assert(v.size()==nlin() && j<ncol());
    #ifdef HAVE_BLAS
        const BLAS_INT M = sizet_to_int(nlin());
        BLAS(dcopy,DCOPY)(M,v.data(),1,data()+nlin()*j,1);
    #else
        for (Index i=0; i<nlin(); ++i)
            (*this)(i,j) = v(i);
    #endif
    }

    inline void Matrix::setlin(const Index i,const Vector& v) {
        om_assert(v.size()==ncol());
        om_assert(i<nlin());
    #ifdef HAVE_BLAS
        const BLAS_INT M = sizet_to_int(nlin());
        const BLAS_INT N = sizet_to_int(ncol());
        BLAS(dcopy,DCOPY)(N,v.data(),1,data()+i,M);
    #else
        for (Index j=0; j<ncol(); ++j)
            (*this)(i,j) = v(j);
    #endif
    }

    inline Vector Matrix::tmult(const Vector& v) const {
        om_assert(nlin()==v.nlin());
        Vector res(ncol());
    #ifdef HAVE_BLAS
        const BLAS_INT M = sizet_to_int(nlin());
        const BLAS_INT N = sizet_to_int(ncol());
        DGEMV(CblasTrans,M,N,1.0,data(),sizet_to_int(nlin()),v.data(),1,0.0,res.data(),1);
    #else
        for (Index i=0; i<ncol(); ++i) {
            res(i) = 0;
            for (Index j=0; j<nlin(); ++j)
                res(i) += (*this)(j,i)*v(j);
        }
    #endif

        return res;
    }

    inline Matrix Matrix::inverse() const {
        om_assert(nlin()==ncol());
    #ifdef HAVE_LAPACK
        Matrix invA(*this,DEEP_COPY);
        // LU
        #if defined(CLAPACK_INTERFACE)
            const BLAS_INT M = sizet_to_int(invA.nlin());
            const BLAS_INT N = sizet_to_int(ncol());
            BLAS_INT* pivots = new BLAS_INT[N];
            DGETRF(M,N,invA.data(),M,pivots);
            DGETRI(N,invA.data(),N,pivots);
            delete[] pivots;
        #else
            int Info = 0;
            BLAS_INT M = sizet_to_int(invA.nlin());
            BLAS_INT N = sizet_to_int(ncol());
            int* pivots = new int[N];
            DGETRF(M,N,invA.data(),M,pivots,Info);
            const Dimension sz = invA.ncol()*64;
            double* work=new double[sz];
            DGETRI(N,invA.data(),N,pivots,work,sz,Info);
            delete[] pivots;
            delete[] work;
            om_assert(Info==0);
        #endif
        return invA;
    #else
        throw OpenMEEG::maths::LinearAlgebraError("Inverse not implemented, requires LAPACK");
    #endif
    }

    inline Matrix Matrix::operator*(const Matrix& B) const {
        om_assert(ncol()==B.nlin());
        Matrix C(nlin(),B.ncol());
    #ifdef HAVE_BLAS
        const BLAS_INT M = sizet_to_int(nlin());
        const BLAS_INT N = sizet_to_int(ncol());
        const BLAS_INT L = sizet_to_int(B.ncol());
        DGEMM(CblasNoTrans,CblasNoTrans,M,L,N,1.0,data(),M,B.data(),N,0.0,C.data(),M);
    #else
        for (Index i=0; i<C.nlin(); ++i)
            for (Index j=0; j<C.ncol(); ++j) {
                C(i,j) = 0.0;
                for (Index k=0; k<ncol(); ++k)
                    C(i,j) += (*this)(i,k)*B(k,j);
            }
    #endif
        return C;
    }

    inline Matrix Matrix::tmult(const Matrix& B) const {
        om_assert(nlin()==B.nlin());
        Matrix C(ncol(),B.ncol());
    #ifdef HAVE_BLAS
        const BLAS_INT M = sizet_to_int(nlin());
        const BLAS_INT N = sizet_to_int(ncol());
        const BLAS_INT L = sizet_to_int(B.ncol());
        DGEMM(CblasTrans,CblasNoTrans,N,L,M,1.0,data(),M,B.data(),M,0.0,C.data(),N);
    #else
        for (Index i=0; i<C.nlin(); ++i)
            for (Index j=0; j<C.ncol(); ++j) {
                C(i,j) = 0.0;
                for (Index k=0; k<nlin(); ++k)
                    C(i,j) += (*this)(k,i)*B(k,j);
            }
    #endif
        return C;
    }

    inline Matrix Matrix::multt(const Matrix& B) const {
        om_assert(ncol()==B.ncol());
        Matrix C(nlin(),B.nlin());
    #ifdef HAVE_BLAS
        const BLAS_INT M = sizet_to_int(nlin());
        const BLAS_INT N = sizet_to_int(ncol());
        const BLAS_INT L = sizet_to_int(B.nlin());
        DGEMM(CblasNoTrans,CblasTrans,M,L,N,1.0,data(),M,B.data(),L,0.0,C.data(),M);
    #else
        for (Index j=0; j<C.ncol(); ++j)
            for (Index i=0; i<C.nlin(); ++i) {
                C(i,j) = 0.0;
                for (Index k=0; k<ncol(); ++k)
                    C(i,j) += (*this)(i,k)*B(j,k);
            }
    #endif
        return C;
    }

    inline Matrix Matrix::tmultt(const Matrix& B) const {
        om_assert(nlin()==B.ncol());
        Matrix C(ncol(),B.nlin());
    #ifdef HAVE_BLAS
        const BLAS_INT M = sizet_to_int(nlin());
        const BLAS_INT N = sizet_to_int(ncol());
        const BLAS_INT L = sizet_to_int(B.nlin());
        DGEMM(CblasTrans,CblasTrans,L,N,M,1.0,data(),M,B.data(),N,0.0,C.data(),L);
    #else
        for (Index i=0; i<C.nlin(); ++i)
            for (Index j=0; j<C.ncol(); ++j) {
                C(i,j) = 0.0;
                for (Index k=0; k<nlin(); ++k)
                    C(i,j) += (*this)(k,i)*B(j,k);
            }
    #endif
        return C;
    }

    inline Matrix Matrix::operator*(const SymMatrix& B) const {
        om_assert(ncol()==B.nlin());
        Matrix C(nlin(),B.ncol());

    // Workaround an MKL bug
    //#ifdef HAVE_BLAS
    #if defined(HAVE_BLAS) && !defined(USE_MKL)
        Matrix D(B);
        const BLAS_INT m = sizet_to_int(nlin());
        const BLAS_INT n = sizet_to_int(B.ncol());
        DSYMM(CblasRight,CblasUpper,m,n,1.0,D.data(),n,data(),m,0.0,C.data(),m);
    #else
        for (Index j=0; j<B.ncol(); ++j)
            for (Index i=0; i<nlin(); ++i) {
                double sum = 0.0;
                for (size_t k=0; k<ncol(); ++k)
                    sum += (*this)(i,k)*B(k,j);
                C(i,j) = sum;
            }
    #endif
        return C;
    }

    inline void Matrix::operator+=(const Matrix& B) {
        om_assert(nlin()==B.nlin());
        om_assert(ncol()==B.ncol());
    #ifdef HAVE_BLAS
        const BLAS_INT sz = sizet_to_int(size());
        BLAS(daxpy,DAXPY)(sz,1.0,B.data(),1,data(),1);
    #else
        const size_t sz = size();
        for (size_t i=0; i<sz; ++i)
            data()[i] += B.data()[i];
    #endif
    }

    inline void Matrix::operator-=(const Matrix& B) {
        om_assert(nlin()==B.nlin());
        om_assert(ncol()==B.ncol());
    #ifdef HAVE_BLAS
        const BLAS_INT sz = sizet_to_int(size());
        BLAS(daxpy,DAXPY)(sz,-1.0,B.data(),1,data(),1);
    #else
        const size_t sz = size();
        for (size_t i=0; i<sz; ++i)
            data()[i] -= B.data()[i];
    #endif
    }

    inline double Matrix::dot(const Matrix& B) const {
        om_assert(nlin()==B.nlin());
        om_assert(ncol()==B.ncol());
    #ifdef HAVE_BLAS
        const BLAS_INT sz = sizet_to_int(size());
        return BLAS(ddot,DDOT)(sz,data(),1,B.data(),1);
    #else
        const  sz = size();
        double s = 0.0;
        for (size_t i=0; i<sz; i++)
            s += data()[i]*B.data()[i];
        return s;
    #endif
    }
}
