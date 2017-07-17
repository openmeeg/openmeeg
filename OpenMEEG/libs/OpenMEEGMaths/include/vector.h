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

#include <OMassert.H>
#include <cstdlib>
#include <string>

#include <OpenMEEGMathsConfig.h>
#include <linop.h>
#include <RC.H>
#include <MathsIO.H>

namespace OpenMEEG {

    class Matrix;
    class SymMatrix;

    class OPENMEEGMATHS_EXPORT Vector: public LinOp {

        utils::RCPtr<LinOpValue> value;

    public:

        Vector(): LinOp(0,1,FULL,1),value() { }

        Vector(const size_t N): LinOp(N,1,FULL,1),value(new LinOpValue(size())) { }
        Vector(const Vector& A,const DeepCopy): LinOp(A.nlin(),1,FULL,1),value(new LinOpValue(A.size(),A.data())) { }

        explicit Vector(Matrix& A);
        explicit Vector(SymMatrix& A);

        void alloc_data() { value = new LinOpValue(size()); }
        void reference_data(const double* array) { value = new LinOpValue(size(),array); }

        size_t size() const { return nlin(); }

        bool empty() const { return value->empty(); }

        double* data() const { return value->data; }

        inline double operator()(const size_t i) const {
            om_assert(i<nlin());
            return value->data[i];
        }

        inline double& operator()(const size_t i) {
            om_assert(i<nlin());
            return value->data[i];
        }

        Vector subvect(size_t istart, size_t isize) const;
        Vector operator+(const Vector& v) const;
        Vector operator-(const Vector& v) const;
        void operator+=(const Vector& v);
        void operator-=(const Vector& v);
        void operator*=(double x);
        void operator/=(double x) { (*this) *= (1.0/x); }
        Vector operator+(double i) const;
        Vector operator-(double i) const;
        Vector operator*(double x) const;
        Vector operator/(double x) const { return (*this)*(1.0/x); }
        double operator*(const Vector& v) const;
        Vector operator*(const Matrix& m) const;

        Vector kmult(const Vector& x) const;
        // Vector conv(const Vector& v) const;
        // Vector conv_trunc(const Vector& v) const;
        Matrix outer_product(const Vector& v) const;

        double norm() const;
        double sum() const;
        double mean() const { return sum()/size(); }

        void set(double x);
        void save(const char *filename) const;
        void load(const char *filename);

        void save(const std::string& s) const { save(s.c_str()); }
        void load(const std::string& s)       { load(s.c_str()); }

        void info() const;

        friend class SymMatrix;
        friend class Matrix;
    };

    OPENMEEGMATHS_EXPORT Vector operator*(const double &d, const Vector &v);

    OPENMEEGMATHS_EXPORT std::ostream& operator<<(std::ostream& f,const Vector &M);
    OPENMEEGMATHS_EXPORT std::istream& operator>>(std::istream& f,Vector &M);

    inline Vector Vector::subvect(size_t istart, size_t isize) const {
        om_assert (istart+isize<=nlin());
        Vector a(isize);
        for (size_t i=0; i<isize; i++)
            a(i) = (*this)(istart+i);
        return a;
    }

    inline Vector Vector::operator+(const Vector& v) const {
        om_assert(nlin()==v.nlin());
        Vector p(*this,DEEP_COPY);
    #ifdef HAVE_BLAS
        BLAS(daxpy,DAXPY)(sizet_to_int(nlin()),1,v.data(),1,p.data(),1);
    #else
        for( size_t i=0; i<nlin(); i++ )
            p.data()[i]=data()[i]+v.data()[i];
    #endif
        return p;
    }

    inline Vector Vector::operator-(const Vector& v) const {
        om_assert(nlin()==v.nlin());
        Vector p(*this,DEEP_COPY);
    #ifdef HAVE_BLAS
        BLAS(daxpy,DAXPY)(sizet_to_int(nlin()),-1,v.data(),1,p.data(),1);
    #else
        for( size_t i=0; i<nlin(); i++ )
            p.data()[i]=data()[i]-v.data()[i];
    #endif
        return p;
    }

    inline void Vector::operator+=(const Vector& v) {
        om_assert(nlin()==v.nlin());
    #ifdef HAVE_BLAS
        BLAS(daxpy,DAXPY)(sizet_to_int(nlin()),1,v.data(),1,data(),1);
    #else
        for( size_t i=0; i<nlin(); i++ )
            data()[i]+=v.data()[i];
    #endif
    }

    inline void Vector::operator-=(const Vector& v) {
        om_assert(nlin()==v.nlin());
    #ifdef HAVE_BLAS
        BLAS(daxpy,DAXPY)(sizet_to_int(nlin()),-1,v.data(),1,data(),1);
    #else
        for( size_t i=0; i<nlin(); i++ )
            data()[i]-=v.data()[i];
    #endif
    }

    inline double Vector::operator*(const Vector& v) const {
        om_assert(nlin()==v.nlin());
    #ifdef HAVE_BLAS
        return BLAS(ddot,DDOT)(sizet_to_int(nlin()),data(),1,v.data(),1);
    #else
        double s=0;
        for( size_t i=0; i<nlin(); i++ )
            s+=data()[i]*v.data()[i];
        return s;
    #endif
    }

    inline Vector Vector::operator*(double x) const {
    #ifdef HAVE_BLAS
        Vector p(*this,DEEP_COPY);
        BLAS(dscal,DSCAL)(sizet_to_int(nlin()),x,p.data(),1);
    #else
        Vector p(nlin());
        for( size_t i=0; i<nlin(); i++ )
            p.data()[i]=x*data()[i];
    #endif
        return p;
    }

    inline void Vector::operator*=(double x) {
    #ifdef HAVE_BLAS
        BLAS(dscal,DSCAL)(sizet_to_int(nlin()),x,data(),1);
    #else
        for( size_t i=0; i<nlin(); i++ )
            data()[i]*=x;
    #endif
    }

    inline double Vector::norm() const
    {
    #ifdef HAVE_BLAS
        return BLAS(dnrm2,DNRM2)(sizet_to_int(nlin()),data(),1);
    #else
        std::cout << "'Vector::norm' not implemented" << std::endl;
        exit(1);
        return 0;
    #endif
    }

    // inline Vector Vector::conv(const Vector& v) const {
    //     if (v.nlin()<nlin()) return v.conv(*this);
    // 
    //     Vector p(nlin()+v.nlin()-1);
    //     p.set(0);
    //     for (size_t i=0; i<v.nlin(); i++) {
    // #ifdef HAVE_BLAS
    //         BLAS(daxpy,DAXPY)(sizet_to_int(nlin()), v(i), data(), 1, p.data()+i, 1);
    // #else
    //         for (size_t j=0;j<nlin();j++)
    //             p(i+j)+=v(i)*data()[j];
    // #endif
    //     }
    //     return p;
    // }
    // 
    // inline Vector Vector::conv_trunc(const Vector& v) const {
    //     Vector p(v.nlin());
    //     p.set(0);
    //     for (size_t i=0; i<v.nlin(); i++)
    //     {
    //         size_t m = std::min(nlin(),v.nlin()-i);
    // #ifdef HAVE_BLAS
    //         BLAS(daxpy,DAXPY)((int)m, v(i), data(), 1, p.data()+i, 1);
    // #else
    //         for (size_t j=0;j<m;j++)
    //             p(i+j)+=v(i)*data()[j];
    // #endif
    //     }
    //     return p;
    // }

    //  Operators.

    OPENMEEGMATHS_EXPORT inline Vector operator*(const double &d, const Vector &v) { return v*d; }
}
