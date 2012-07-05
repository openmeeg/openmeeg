/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
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

#ifndef OPENMEEG_DIAG_MATRIX_H
#define OPENMEEG_DIAG_MATRIX_H

#include <cassert>
#include <fstream>
#include <iostream>
#include <utility>
#include <string>

#include <om_utils.h>
#include <linop.h>
#include <vector.h>

namespace OpenMEEG {

    class OPENMEEGMATHS_EXPORT DiagMatrix : public Vector {

    public:
        DiagMatrix(): Vector() {}

        // DiagMatrix(const char* fname): Vector(fname) {}
        DiagMatrix(size_t N):          Vector(N) {}
        DiagMatrix(size_t M,size_t N): Vector(M) {assert(M==N); }
        DiagMatrix(const DiagMatrix& D,const DeepCopy): Vector(D,DEEP_COPY) { }

        ~DiagMatrix() {};

        inline double operator()( size_t i, size_t j ) const {
            assert(i < nlin() && j < ncol());
            if(i==j)
                return data()[i];
            else
                return 0.;
        }

        inline double& operator()( size_t i, size_t j ) {
            assert(i < nlin() && i == j);
            return data()[i];
        }

        inline double operator()( size_t i ) const { return (*this)(i,i);}
        inline double& operator()( size_t i ) { return (*this)(i,i);}

        inline Vector operator*(const Vector& v) const ;

        size_t ncol() const {return nlin();}
        inline void info() const;

        DiagMatrix inverse() const {
            DiagMatrix inv(nlin());
            for (size_t i=0;i<size();i++)
                inv.data()[i]=1./this->data()[i];
            return inv;
        }
    };

    inline Vector DiagMatrix::operator*(const Vector& v) const {
        Vector result(nlin());
        for (size_t i=0;i<size();i++)
            result(i)=v(i)*this->data()[i];
        return result;
    }

    inline void DiagMatrix::info() const {
        if (size() == 0) {
            std::cout << "DiagMatrix Empty" << std::endl;
            return;
        }

        std::cout << "Size : " << nlin() << std::endl;

        double minv = this->operator()(0,0);
        double maxv = this->operator()(0,0);
        size_t mini = 0;
        size_t maxi = 0;

        for(size_t i = 0; i < nlin(); ++i)
            if (minv > this->operator()(i,i)) {
                minv = this->operator()(i,i);
                mini = i;
            } else if (maxv < this->operator()(i,i)) {
                maxv = this->operator()(i,i);
                maxi = i;
            }

        std::cout << "Min Value : " << minv << " (" << mini << ")" << std::endl;
        std::cout << "Max Value : " << maxv << " (" << maxi << ")" << std::endl;
        std::cout << "First Values" << std::endl;
        for(size_t i = 0; i < std::min(nlin(),(size_t) 5); ++i) {
            for(size_t j = 0; j < std::min(nlin(),(size_t) 5); ++j) 
            {
                std::cout << this->operator()(i,j) << " ";
            }
            std::cout << std::endl;
        }
    }
}

#endif //! OPENMEEG_DIAG_MATRIX_H

