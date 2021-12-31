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
            const double value = (*this)(i);
            if (minv>value) {
                minv = value;
                mini = i;
            } else if (maxv<value ) {
                maxv = value;
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
        catch (maths::Exception& e) {
            ifs >> *this;
        }
    }

    void Vector::save(const char* filename) const {
        maths::ofstream ofs(filename);
        try {
            ofs << maths::format(filename,maths::format::FromSuffix) << *this;
        }
        catch (maths::Exception& e) {
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
