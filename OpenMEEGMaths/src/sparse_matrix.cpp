// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include "sparse_matrix.h"
#include "symmatrix.h"

namespace OpenMEEG {

    double SparseMatrix::frobenius_norm() const {
        double d = 0.;
        for ( const_iterator it = m_tank.begin() ; it != m_tank.end(); ++it) {
            d += std::pow(it->second,2);
        }
        return sqrt(d);
    }

    Vector SparseMatrix::operator*(const Vector &x) const
    {
        Vector ret(nlin());
        ret.set(0);

        Tank::const_iterator it;
        for(it = m_tank.begin(); it != m_tank.end(); ++it) {
            size_t i = it->first.first;
            size_t j = it->first.second;
            double val = it->second;
            ret(i) += val * x(j);
        }

        return ret;
    }

    Matrix SparseMatrix::operator*(const SymMatrix &mat) const
    {
        om_assert(ncol()==mat.nlin());
        Matrix out(nlin(),mat.ncol());
        out.set(0.0);

        Tank::const_iterator it;
        for(it = m_tank.begin(); it != m_tank.end(); ++it) {
            size_t i = it->first.first;
            size_t j = it->first.second;
            double val = it->second;
            for(size_t k = 0; k < mat.ncol(); ++k) {
                out(i,k) += val * mat(j,k);
            }
        }

        return out;
    }

    Matrix SparseMatrix::operator*(const Matrix &mat) const
    {
        om_assert(ncol()==mat.nlin());
        Matrix out(nlin(),mat.ncol());
        out.set(0.0);

        for( Tank::const_iterator it = m_tank.begin(); it != m_tank.end(); ++it) {
            size_t i = it->first.first;
            size_t j = it->first.second;
            double val = it->second;
            for ( size_t k = 0; k < mat.ncol(); ++k) {
                out(i, k) += val * mat(j, k);
            }
        }

        return out;
    }

    SparseMatrix SparseMatrix::operator*(const SparseMatrix &mat) const
    {
        // fast enough ?
        om_assert(ncol() == mat.nlin());
        SparseMatrix out(nlin(), mat.ncol());

        for ( Tank::const_iterator it1 = m_tank.begin(); it1 != m_tank.end(); ++it1) {
            size_t i = it1->first.first;
            size_t j = it1->first.second;
            for ( Tank::const_iterator it2 = mat.begin(); it2 != mat.end(); ++it2) {
                if ( it2->first.first == j ) {
                    out(i, it2->first.second) += it1->second * it2->second;
                }
            }
        }
        return out;
    }

    SparseMatrix SparseMatrix::operator+(const SparseMatrix &mat) const
    {
        om_assert(nlin() == mat.nlin() && ncol() == mat.ncol());
        SparseMatrix out(nlin(), ncol());

        for ( Tank::const_iterator it = m_tank.begin(); it != m_tank.end(); ++it) {
            size_t i = it->first.first;
            size_t j = it->first.second;
            out(i, j) += it->second;
        }
        for ( Tank::const_iterator it = mat.begin(); it != mat.end(); ++it) {
            size_t i = it->first.first;
            size_t j = it->first.second;
            out(i, j) += it->second;
        }
        return out;
    }

    SparseMatrix SparseMatrix::transpose() const {
        SparseMatrix tsp(ncol(),nlin());
        const_iterator it;
        for(it = m_tank.begin(); it != m_tank.end(); ++it) {
            size_t i = it->first.first;
            size_t j = it->first.second;
            tsp(j,i) = it->second;
        }
        return tsp;
    }

    void SparseMatrix::set(double d) {
        for( iterator it = m_tank.begin(); it != m_tank.end(); ++it) {
            it->second = d;
        }
    }

    void SparseMatrix::info() const {
        if ((nlin() == 0) || (ncol() == 0) || m_tank.empty()) {
            std::cout << "Matrix Empty" << std::endl;
            return;
        }

        std::cout << "Dimensions : " << nlin() << " x " << ncol() << std::endl;

        double minv = m_tank.begin()->second;
        double maxv = m_tank.begin()->second;
        size_t mini = 0;
        size_t maxi = 0;
        size_t minj = 0;
        size_t maxj = 0;

        for (Tank::const_iterator it=m_tank.begin();it!=m_tank.end();++it) {
            if (minv>it->second) {
                minv = it->second;
                mini = it->first.first;
                minj = it->first.second;
            } else if (maxv<it->second) {
                maxv = it->second;
                maxi = it->first.first;
                maxj = it->first.second;
            }
        }

        std::cout << "Min Value : " << minv << " (" << mini << "," << minj << ")" << std::endl;
        std::cout << "Max Value : " << maxv << " (" << maxi << "," << maxj << ")" << std::endl;
        std::cout << "First Values" << std::endl;

        size_t cnt = 0;
        for(Tank::const_iterator it = m_tank.begin(); it != m_tank.end() && cnt < 5; ++it) {
            std::cout << "(" << it->first.first << "," << it->first.second << ") " << it->second << std::endl;
            cnt++;
        }
    }

    // =======
    // = IOs =
    // =======

    void SparseMatrix::load(const char *filename) {
        maths::ifstream ifs(filename);
        try {
            ifs >> maths::format(filename,maths::format::FromSuffix) >> *this;
        }
        catch (maths::Exception& e) {
            ifs >> *this;
        }
    }

    void SparseMatrix::save(const char *filename) const {
        maths::ofstream ofs(filename);
        try {
            ofs << maths::format(filename,maths::format::FromSuffix) << *this;
        }
        catch (maths::Exception& e) {
            ofs << *this;
        }
    }
}
