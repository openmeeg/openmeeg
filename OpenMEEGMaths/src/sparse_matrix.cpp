// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include "sparse_matrix.h"
#include "symmatrix.h"

namespace OpenMEEG {

    static inline double sqr(const double x) { return x*x; }

    double SparseMatrix::frobenius_norm() const {
        double d = 0.;
        for (const auto& tkelmt : m_tank)
            d += sqr(tkelmt.second);
        return sqrt(d);
    }

    Vector SparseMatrix::operator*(const Vector& x) const {
        Vector ret(nlin());
        ret.set(0);

        for (const auto& tkelmt : m_tank) {
            const size_t i = tkelmt.first.first;
            const size_t j = tkelmt.first.second;
            const double val = tkelmt.second;
            ret(i) += val*x(j);
        }

        return ret;
    }

    Matrix SparseMatrix::operator*(const SymMatrix& mat) const
    {
        om_assert(ncol()==mat.nlin());
        Matrix out(nlin(),mat.ncol());
        out.set(0.0);

        for (const auto& tkelmt : m_tank) {
            const size_t i = tkelmt.first.first;
            const size_t j = tkelmt.first.second;
            const double val = tkelmt.second;
            for (size_t k=0; k<mat.ncol(); ++k)
                out(i,k) += val*mat(j,k);
        }

        return out;
    }

    Matrix SparseMatrix::operator*(const Matrix& mat) const {
        om_assert(ncol()==mat.nlin());
        Matrix out(nlin(),mat.ncol());
        out.set(0.0);

        for (const auto& tkelmt : m_tank) {
            const size_t i = tkelmt.first.first;
            const size_t j = tkelmt.first.second;
            const double val = tkelmt.second;
            for (size_t k=0; k<mat.ncol(); ++k)
                out(i,k) += val*mat(j,k);
        }

        return out;
    }

    SparseMatrix SparseMatrix::operator*(const SparseMatrix& mat) const {
        // fast enough ?
        om_assert(ncol()==mat.nlin());
        SparseMatrix out(nlin(),mat.ncol());

        for (const auto& tkelmt1 : m_tank) {
            const size_t i = tkelmt1.first.first;
            const size_t j = tkelmt1.first.second;
            for (const auto& elmt2 : mat)
                if (elmt2.first.first==j)
                    out(i,elmt2.first.second) += tkelmt1.second*elmt2.second;
        }
        return out;
    }

    SparseMatrix SparseMatrix::operator+(const SparseMatrix& mat) const {
        om_assert(nlin()==mat.nlin() && ncol()==mat.ncol());
        SparseMatrix out(nlin(),ncol());

        for (const auto& tkelmt : m_tank) {
            const size_t i = tkelmt.first.first;
            const size_t j = tkelmt.first.second;
            out(i,j) += tkelmt.second;
        }
        for (const auto& tkelmt : m_tank) {
            const size_t i = tkelmt.first.first;
            const size_t j = tkelmt.first.second;
            out(i,j) += tkelmt.second;
        }
        return out;
    }

    SparseMatrix SparseMatrix::transpose() const {
        SparseMatrix tsp(ncol(),nlin());
        for (const auto& tkelmt : m_tank) {
            const size_t i = tkelmt.first.first;
            const size_t j = tkelmt.first.second;
            tsp(j,i) = tkelmt.second;
        }
        return tsp;
    }


    void SparseMatrix::info() const {
        if (nlin()==0 || ncol()==0 || m_tank.empty()) {
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

        for (const auto& tkelmt : m_tank)
            if (minv>tkelmt.second) {
                minv = tkelmt.second;
                mini = tkelmt.first.first;
                minj = tkelmt.first.second;
            } else if (maxv<tkelmt.second) {
                maxv = tkelmt.second;
                maxi = tkelmt.first.first;
                maxj = tkelmt.first.second;
            }

        std::cout << "Min Value : " << minv << " (" << mini << "," << minj << ")" << std::endl;
        std::cout << "Max Value : " << maxv << " (" << maxi << "," << maxj << ")" << std::endl;
        std::cout << "First Values" << std::endl;

        size_t cnt = 0;
        for (const auto& tkelmt : m_tank) {
            std::cout << "(" << tkelmt.first.first << "," << tkelmt.first.second << ") " << tkelmt.second << std::endl;
            if (++cnt==5)
                break;
        }
    }
}
