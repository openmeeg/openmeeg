// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <OMassert.H>
#include <cstdlib>
#include <string>

#include <OpenMEEGMathsConfig.h>
#include <linop.h>
#include <MathsIO.H>

namespace OpenMEEG {

    class Matrix;
    class SymMatrix;

    class OPENMEEGMATHS_EXPORT Vector: public LinOp {

        LinOpValue value;

    public:

        Vector(): LinOp(0,1,FULL,1),value() { }

        Vector(const Dimension N): LinOp(N,1,FULL,1),value(N) { }
        Vector(const Vector& A,const DeepCopy): LinOp(A.nlin(),1,FULL,1),value(A.size(),A.data()) { }

        explicit Vector(const Matrix& A);
        explicit Vector(const SymMatrix& A);

        void alloc_data() { value = LinOpValue(size()); }
        void reference_data(const double* array) { value = LinOpValue(size(),array); }

        size_t size() const { return nlin(); }

        bool empty() const { return value.empty(); }

        double* data() const { return value.get(); }

        double operator()(const Index i) const {
            om_assert(i<nlin());
            return value[i];
        }

        double& operator()(const Index i) {
            om_assert(i<nlin());
            return value[i];
        }

        Vector subvect(const Index istart,const Index isize) const;

        Vector operator+(const Vector& v) const;
        Vector operator-(const Vector& v) const;

        Vector operator-() const {
            // No Blas implementation ?
            Vector res(nlin());
            for (Index i=0; i<nlin(); i++ )
                res.data()[i] = -data()[i];
            return res;
        }

        inline void operator+=(const Vector& v);
        inline void operator-=(const Vector& v);
        inline void operator*=(const double x);
        void operator/=(const double x) { (*this) *= (1.0/x); }
        Vector operator+(const double i) const;
        Vector operator-(const double i) const;
        inline Vector operator*(const double x) const;
        Vector operator/(const double x) const { return (*this)*(1.0/x); }
        inline double operator*(const Vector& v) const;
        Vector operator*(const Matrix& m) const;

        Vector kmult(const Vector& x) const;
        // Vector conv(const Vector& v) const;
        // Vector conv_trunc(const Vector& v) const;
        Matrix outer_product(const Vector& v) const;

        double norm() const;
        double sum() const;
        double mean() const { return sum()/size(); }

        void set(const double x);

        void save(const char *filename) const;
        void load(const char *filename);

        void save(const std::string& s) const { save(s.c_str()); }
        void load(const std::string& s)       { load(s.c_str()); }

        void info() const;

        friend class SymMatrix;
        friend class Matrix;
    };

    OPENMEEGMATHS_EXPORT Vector operator*(const double d,const Vector& v);

    OPENMEEGMATHS_EXPORT std::ostream& operator<<(std::ostream& f,const Vector& M);
    OPENMEEGMATHS_EXPORT std::istream& operator>>(std::istream& f,Vector& M);

    inline Vector Vector::subvect(const Index istart,const Index isize) const {
        om_assert (istart+isize<=nlin());
        Vector a(isize);
        for (Index i=0; i<isize; ++i)
            a(i) = (*this)(istart+i);
        return a;
    }

    inline Vector Vector::operator+(const Vector& v) const {
        om_assert(nlin()==v.nlin());
        Vector p(*this,DEEP_COPY);
    #ifdef HAVE_BLAS
        BLAS(daxpy,DAXPY)(sizet_to_int(nlin()),1,v.data(),1,p.data(),1);
    #else
        for (Index i=0; i<nlin(); ++i)
            p.data()[i] = data()[i]+v.data()[i];
    #endif
        return p;
    }

    inline Vector Vector::operator-(const Vector& v) const {
        om_assert(nlin()==v.nlin());
        Vector p(*this,DEEP_COPY);
    #ifdef HAVE_BLAS
        BLAS(daxpy,DAXPY)(sizet_to_int(nlin()),-1,v.data(),1,p.data(),1);
    #else
        for (Index i=0; i<nlin(); ++i)
            p.data()[i] = data()[i]-v.data()[i];
    #endif
        return p;
    }

    inline void Vector::operator+=(const Vector& v) {
        om_assert(nlin()==v.nlin());
    #ifdef HAVE_BLAS
        BLAS(daxpy,DAXPY)(sizet_to_int(nlin()),1,v.data(),1,data(),1);
    #else
        for (Index i=0; i<nlin(); ++i)
            data()[i] += v.data()[i];
    #endif
    }

    inline void Vector::operator-=(const Vector& v) {
        om_assert(nlin()==v.nlin());
    #ifdef HAVE_BLAS
        BLAS(daxpy,DAXPY)(sizet_to_int(nlin()),-1,v.data(),1,data(),1);
    #else
        for (Index i=0; i<nlin(); ++i)
            data()[i] -= v.data()[i];
    #endif
    }

    inline double Vector::operator*(const Vector& v) const {
        om_assert(nlin()==v.nlin());
    #ifdef HAVE_BLAS
        return BLAS(ddot,DDOT)(sizet_to_int(nlin()),data(),1,v.data(),1);
    #else
        double s=0;
        for (Index i=0; i<nlin(); ++i)
            s += data()[i]*v.data()[i];
        return s;
    #endif
    }

    inline Vector Vector::operator*(const double x) const {
    #ifdef HAVE_BLAS
        Vector p(*this,DEEP_COPY);
        BLAS(dscal,DSCAL)(sizet_to_int(nlin()),x,p.data(),1);
    #else
        Vector p(nlin());
        for (Index i=0; i<nlin(); ++i)
            p.data()[i] = x*data()[i];
    #endif
        return p;
    }

    inline void Vector::operator*=(const double x) {
    #ifdef HAVE_BLAS
        BLAS(dscal,DSCAL)(sizet_to_int(nlin()),x,data(),1);
    #else
        for (Index i=0; i<nlin(); ++i)
            data()[i] *= x;
    #endif
    }

    inline double Vector::norm() const {
    #ifdef HAVE_BLAS
        return BLAS(dnrm2,DNRM2)(sizet_to_int(nlin()),data(),1);
    #else
        throw OpenMEEG::maths::LinearAlgebraError("'Vector::norm' not implemented, requires BLAS");
    #endif
    }

    // inline Vector Vector::conv(const Vector& v) const {
    //     if (v.nlin()<nlin())
    //         return v.conv(*this);
    //
    //     Vector p(nlin()+v.nlin()-1);
    //     p.set(0);
    //     for (Index i=0; i<v.nlin(); ++i) {
    // #ifdef HAVE_BLAS
    //         BLAS(daxpy,DAXPY)(sizet_to_int(nlin()),v(i),data(),1,p.data()+i,1);
    // #else
    //         for (Index j=0; j<nlin(); ++j)
    //             p(i+j) += v(i)*data()[j];
    // #endif
    //     }
    //     return p;
    // }
    //
    // inline Vector Vector::conv_trunc(const Vector& v) const {
    //     Vector p(v.nlin());
    //     p.set(0);
    //     for (Index i=0; i<v.nlin(); ++i)
    //     {
    //         const Index m = std::min(nlin(),v.nlin()-i);
    // #ifdef HAVE_BLAS
    //         BLAS(daxpy,DAXPY)((int)m,v(i),data(),1,p.data()+i,1);
    // #else
    //         for (Index j=0; j<m; ++j)
    //             p(i+j) += v(i)*(*this)(j);
    // #endif
    //     }
    //     return p;
    // }

    //  Operators.

    inline Vector operator*(const double d,const Vector& v) { return v*d; }
}
