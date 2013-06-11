#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <cmath>
#include <lapack++.H>

namespace Maths {
    
    template <typename TYPE>
    class Vector {
    public:

        Vector(const unsigned d): n(d),elmts(new TYPE[d]) { }
        ~Vector() { delete[] elmts; }

        Vector(const Vector& V): n(V.n),elmts(new TYPE[n]) { std::copy(V.elmts,V.elmts+n,elmts); }

              TYPE& operator()(const unsigned int i)      { return elmts[i-1]; }
        const TYPE&  operator()(const unsigned int i) const{ return elmts[i-1]; }

        unsigned dim() const { return n; }

        Vector& operator=(const Vector<TYPE>& V) {
            for (unsigned i=1;i<=dim();++i)
                (*this)(i) = V(i);
            return *this;
        }

        Vector& operator=(const TYPE val) {
            for (unsigned i=1;i<=dim();++i)
                (*this)(i) = val;
            return *this;
        }

        Vector& operator+=(const Vector<TYPE>& V) {
            for (unsigned i=1;i<=dim();++i)
                (*this)(i) += V(i);
            return *this;
        }

        Vector operator+(const Vector<TYPE>& V) const {
            Vector res = *this;
            return res += V;
        }

        Vector operator-() const {
            Vector res = *this;
            for (unsigned i=1;i<=dim();++i)
                res(i) = -(*this)(i);
            return res;
        }

        Vector& operator-=(const Vector<TYPE>& V) {
            for (unsigned i=1;i<=dim();++i)
                (*this)(i) -= V(i);
            return *this;
        }

        Vector& operator*=(const TYPE val) {
            for (unsigned i=1;i<=dim();++i)
                (*this)(i) *= val;
            return *this;
        }

        Vector& operator/=(const TYPE val) {
            *this *= 1.0/val;
            return *this;
        }

        Vector operator*(const TYPE val) const {
            Vector res = *this;
            return res *= val;
        }

        Vector operator/(const TYPE val) const {
            Vector res = *this;
            return res /= val;
        }

    private:

        unsigned n;
        TYPE*    elmts;
    };

    template <typename TYPE>
    inline Vector<TYPE> operator*(const TYPE s,const Vector<TYPE>& V) {
        return V.operator*(s);
    }

    template <typename TYPE>
    inline TYPE dotprod(const Vector<TYPE>& V1,const Vector<TYPE>& V2) {
        assert(V1.dim()==V2.dim());
        TYPE dp = 0;
        for (unsigned i=1;i<=V1.dim();++i)
            dp += V1(i)*V2(i);
        return dp;
    }

    template <typename TYPE>
    inline TYPE norm(const Vector<TYPE>& V) {
        return std::sqrt(dotprod(V,V));
    }

    template <typename TYPE>
    std::istream& operator>>(std::istream& is,Vector<TYPE>& V) {
        for (unsigned i=1;i<=V.dim();++i)
            is >> V(i);
        return is;
    }

    template <typename TYPE>
    std::ostream& operator<<(std::ostream& os,const Vector<TYPE>& V) {
        for (unsigned i=1;i<=V.dim();++i)
            os << V(i) << ' ';
        return os;
    }

    template <typename TYPE>
    class Matrix {

        void alloc(const unsigned n,const unsigned m) {
            if (n*m<nl*nc)
                return;
            if (elmts!=0)
                delete[]  elmts;
            elmts = new TYPE[n*m];
        }

    public:

        Matrix(): nl(0),nc(0),elmts(0) { }

        Matrix(const unsigned n,const unsigned m): nl(n),nc(m),elmts(new TYPE[n*m]) { }
        ~Matrix() { delete[] elmts; }

        Matrix(const Matrix& M): nl(M.nl),nc(M.nc),elmts(new TYPE[nl*nc]) { std::copy(M.elmts,M.elmts+nl*nc,elmts); }

              TYPE& operator()(const unsigned int i,const unsigned j)      { return elmts[(j-1)*nl+(i-1)]; }
        const TYPE&  operator()(const unsigned int i,const unsigned j) const{ return elmts[(j-1)*nl+(i-1)]; }

        Matrix& operator=(const Matrix& M) {
            if (&M==this)
                return *this;
            alloc(M.rows(),M.cols());
            nl = M.rows();
            nc = M.cols();
            std::copy(M.elmts,M.elmts+nl*nc,elmts);
        }

        Matrix& operator=(const TYPE val) {
            for (unsigned i=1;i<=rows();++i)
                for (unsigned j=1;j<=cols();++j)
                    (*this)(i,j) = val;
            return *this;
        }

        Matrix operator-() const {
            Matrix res = *this;
            for (unsigned i=1;i<=res.rows();++i)
                for (unsigned j=1;j<=res.cols();++j)
                    res(i,j) = - (*this)(i,j);
            return res;
        }

        Vector<TYPE> operator*(const Vector<TYPE>& V) const {
            assert(V.dim()==cols());
            Vector<TYPE> res(rows());
            for (unsigned i=1;i<=rows();++i) {
                TYPE val =0.0;
                for (unsigned j=1;j<=cols();++j)
                    val += (*this)(i,j)*V(j);
                res(i) = val;
            }
            return res;
        }

        Matrix<TYPE> operator*(const Matrix<TYPE>& M) const {
            assert(M.rows()==cols());
            Matrix<TYPE> res(rows(),M.cols());
            for (unsigned i=1;i<=res.rows();++i)
                for (unsigned j=1;j<=res.cols();++j) {
                    TYPE val =0.0;
                    for (unsigned k=1;k<=cols();++k)
                        val += (*this)(i,k)*M(k,j);
                    res(i,j) = val;
                }
            return res;
        }

        Matrix<TYPE>& operator*=(const TYPE val) {
            for (unsigned i=1;i<=rows();++i)
                for (unsigned j=1;j<=cols();++j)
                    (*this)(i,j) *= val;
            return *this;
        }

        Matrix<TYPE>& operator/=(const TYPE val) {
            return operator*=(1.0/val);
        }

        Vector<TYPE> row(const unsigned i) const {
            Vector<TYPE> res(cols());
            for (unsigned j=1;j<=cols();++j)
                res(j) = (*this)(i,j);
            return res;
        }

        Vector<TYPE> col(const unsigned j) const {
            Vector<TYPE> res(rows());
            for (unsigned i=1;i<=rows();++i)
                res(i) = (*this)(i,j);
            return res;
        }

        unsigned rows() const { return nl; }
        unsigned cols() const { return nc; }

    private:

        unsigned nl;
        unsigned nc;
        TYPE*    elmts;
    };

    template <typename TYPE>
    std::ostream& operator<<(std::ostream& os,const Matrix<TYPE>& M) {
        for (unsigned i=1;i<=M.rows();++i) {
            for (unsigned j=1;j<=M.cols();++j)
                os << M(i,j) << ' ';
            os << std::endl;
        }
        return os;
    }

    template <typename TYPE>
    std::istream& operator>>(std::istream& is,Matrix<TYPE>& M) {
        for (unsigned i=1;i<=M.rows();++i) {
            for (unsigned j=1;j<=M.cols();++j)
                is >> M(i,j);
        }
        return is;
    }

    template <typename TYPE>
    inline Matrix<TYPE> transpose(const Matrix<TYPE>& M) {
        Matrix<TYPE> res(M.cols(),M.rows());
        for (unsigned i=1;i<=M.rows();++i)
            for (unsigned j=1;j<=M.cols();++j)
                res(j,i) = M(i,j);
        return res;
    }

    template <typename TYPE> inline
    Matrix<TYPE> inverse(const Matrix<TYPE>& A) {
        assert(A.rows()==A.cols());
        Matrix<TYPE> M = A;
        Lapack<TYPE>::inverse(M.rows(),&M(1,1));
        return M;
    }

    //  Svd Decomposition of matrix A (MxN): A = U*D*V^t.
    //  The singular values of the diagonal matrix D are returned by the function.
    //  These singular values are sorted from largest to smallest.
    //  ### Works only for M>=N ###
    //  If M>N, dimU=[MxN], dimV=NxN, dimD=NxN.

    template <typename TYPE> inline
    Vector<TYPE> Svd(const Matrix<TYPE>& A,Matrix<TYPE>& U,Matrix<TYPE>& V) {
        const unsigned M = A.rows();
        const unsigned N = A.cols();
        Vector<TYPE> D(N);
        Matrix<TYPE> UT = transpose(A);
        Lapack<TYPE>::FullSvd(M,N,&UT(1,1),&D(1),&V(1,1));
        U = transpose(UT);
        return D;
    }

    //  Cholesky decomposition.

    template <typename TYPE> inline
    Matrix<TYPE> Cholesky(const Matrix<TYPE>& A) {
        const unsigned M = A.rows();
        const unsigned N = A.cols();
        assert(M==N);
        Matrix<TYPE> mat(A);
        Lapack<TYPE>::cholesky(N,&mat(1,1));
        return mat;
    }

    //  RQ decomposition.

    template <typename TYPE> inline
    void RQ(const Matrix<TYPE>& A,Matrix<TYPE>& Q,Matrix<TYPE>& R) {
        const unsigned M = A.rows();
        const unsigned N = A.cols();
        Q = A;
        R = 0.0;
        Lapack<TYPE>::RQ(M,N,&Q(1,1),&R(1,1));
    }

    template <typename TYPE> inline
    TYPE det(const Matrix<TYPE>& A) {
        assert(A.rows()==A.cols());
        TYPE MM[A.rows()*A.cols()];
        memcpy(MM,&(A(1,1)),sizeof(MM));
        return Lapack<TYPE>::det(A.rows(),MM);
    }
}

#endif  //  ! MATRIX_H
