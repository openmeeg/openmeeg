// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <OpenMEEGMathsConfig.h>
#include <vector.h>
#include <matrix.h>
#include <symmatrix.h>
#include <algorithm>

namespace OpenMEEG {

    Vector::Vector(const Matrix& A) {
        nlin() = A.nlin()*A.ncol();
        value  = A.value;
    }

    Vector::Vector(const SymMatrix& A) {
        nlin() = A.nlin()*(A.nlin()+1)/2;
        value  = A.value;
    }

    Vector Vector::kmult(const Vector& v) const { // Kronecker multiplication
        om_assert(nlin()==v.nlin());
        Vector p(nlin());
        for (Index i=0; i<nlin(); i++ )
            p(i) = v(i)*(*this)(i);
        return p;
    }

    Vector Vector::operator+(const double x) const {
        Vector p(*this,DEEP_COPY);
        for (Index i=0; i<nlin(); ++i)
            p(i) += x;
        return p;
    }

    Vector Vector::operator-(double x) const {
        Vector p(*this,DEEP_COPY);
        for (Index i=0; i<nlin(); ++i)
            p(i) -= x;
        return p;
    }

    Vector Vector::operator*(const Matrix& m) const {
        om_assert(nlin()==m.nlin());
        Vector c(m.ncol());
        return m.transpose()*(*this);
    }

    void Vector::set(const double x) {
        om_assert(nlin()>0);
        for (Index i=0; i<nlin(); ++i)
            (*this)(i) = x;
    }

    double Vector::sum() const {
        double s=0;
        for (Index i=0; i<nlin(); ++i)
            s += (*this)(i);
        return s;
    }

    void Vector::info() const {
        if (size()==0) {
            std::cout << "Vector Empty" << std::endl;
            return;
        }

        std::cout << "Size : " << nlin() << std::endl;

        double minv = (*this)(0);
        double maxv = (*this)(0);
        Index mini = 0;
        Index maxi = 0;

        for (Index i=0; i<nlin(); ++i) {
            const double val = (*this)(i);
            if (minv>val) {
                minv = val;
                mini = i;
            } else if (maxv<val ) {
                maxv = val;
                maxi = i;
            }
        }

        std::cout << "Min Value : " << minv << " (" << mini << ")" << std::endl;
        std::cout << "Max Value : " << maxv << " (" << maxi << ")" << std::endl;
        std::cout << "First Values" << std::endl;
        for (Index i=0; i<std::min(nlin(),5U); ++i) {
            std::cout << (*this)(i) << std::endl;
        }
    }

    // =======
    // = IOs =
    // =======

    std::ostream& operator<<(std::ostream& f,const Vector& M) {
        for (size_t i=0; i<M.size(); ++i)
            f << M(i) << ' ';
        return f;
    }

    std::istream& operator>>(std::istream& f,Vector &M) {
        for (size_t i=0; i<M.size(); ++i)
            f >> M(i);
        return f;
    }

    void Vector::load(const char* filename) {
        maths::ifstream ifs(filename);
        try {
            ifs >> maths::format(filename, maths::format::FromSuffix) >> *this;
        }
        catch (maths::Exception&) {
            ifs >> *this;
        }
    }

    void Vector::save(const char* filename) const {
        maths::ofstream ofs(filename);
        try {
            ofs << maths::format(filename,maths::format::FromSuffix) << *this;
        }
        catch (maths::Exception&) {
            ofs << *this;
        }
    }

    Matrix Vector::outer_product(const Vector& v) const
    {
        om_assert(size()==v.size());
        Matrix A(size(),v.size());
        A.set(0.0);
    #ifdef HAVE_BLAS
        const BLAS_INT sz = sizet_to_int(size());
        DGER(sz,sz,1.0,data(),1,v.data(),1,A.data(),sz);
    #else
        for(Index j=0; j<nlin(); ++j)
            for(Index i=0; i<nlin(); ++i)
                A(i,j) = v(i)*(*this)(j);
    #endif
        return A;
    }
}
