/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Maureen.Clerc.AT.sophia.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.sophia.inria.fr)

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

namespace OpenMEEG {

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

    Matrix SparseMatrix::operator*(const Matrix &mat) const
    {
        assert(ncol()==mat.nlin());
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

        Tank::const_iterator it;
        for(it = m_tank.begin(); it != m_tank.end(); ++it) {
                if (minv > it->second) {
                    minv = it->second;
                    mini = it->first.first;
                    minj = it->first.second;
                } else if (maxv < it->second) {
                    maxv = it->second;
                    maxi = it->first.first;
                    maxj = it->first.second;
                }
        }

        std::cout << "Min Value : " << minv << " (" << mini << "," << minj << ")" << std::endl;
        std::cout << "Max Value : " << maxv << " (" << maxi << "," << maxj << ")" << std::endl;
        std::cout << "First Values" << std::endl;

        size_t cnt = 0;
        for(it = m_tank.begin(); it != m_tank.end() && cnt < 5; ++it) {
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
            ifs >> maths::format(filename,maths::format::FromSuffix) << *this;
        }
        catch (maths::Exception& e) {
            std::cout << e.what() << " Doing my best...." << std::endl;
            ifs >> *this;
        }
    }

    void SparseMatrix::save(const char *filename) const {
        maths::ofstream ofs(filename);
        try {
            ofs << maths::format(filename,maths::format::FromSuffix) << *this;
        }
        catch (maths::Exception& e) {
            std::cout << e.what() << " Doing my best...." << std::endl;
            ofs << *this;
        }
    }
}
