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
