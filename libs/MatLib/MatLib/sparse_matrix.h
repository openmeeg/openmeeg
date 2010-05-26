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

#ifndef OPENMEEG_SPARSE_MATRIX_H
#define OPENMEEG_SPARSE_MATRIX_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <fstream>
#include <iostream>
#include <map>
#include <utility>

#include "om_utils.h"
#include "linop.h"
#include "vector.h"
#include "matrix.h"

// #ifdef WIN32
// #pragma warning( disable : 4251)    //MSVC warning C4251 : DLL exports of STL templates
// #endif

//template class OPENMEEGMATHS_EXPORT std::pair< size_t, size_t >;
#ifdef WIN32
    template class OPENMEEGMATHS_EXPORT std::map< std::pair< size_t, size_t >, double >;
#endif

namespace OpenMEEG {

    class OPENMEEGMATHS_EXPORT SparseMatrix : public LinOp {

    public:

        typedef std::map< std::pair< size_t, size_t >, double > Tank;
        typedef std::map< std::pair< size_t, size_t >, double >::const_iterator const_iterator;
        typedef std::map< std::pair< size_t, size_t >, double >::iterator iterator;

        SparseMatrix() : LinOp(0,0,SPARSE,2) {};
        SparseMatrix(size_t N,size_t M) : LinOp(N,M,SPARSE,2) {};
        ~SparseMatrix() {};

        inline double operator()( size_t i, size_t j ) const {
            assert(i < nlin());
            assert(j < ncol());
            const_iterator it = m_tank.find(std::make_pair<size_t, size_t>(i, j));
            if (it != m_tank.end()) return it->second;
            else return 0.0;
        }

        inline double& operator()( size_t i, size_t j ) {
            assert(i < nlin());
            assert(j < ncol());
            return m_tank[ std::make_pair( i, j ) ];
        }

        size_t size() const {
            return m_tank.size();
        }

        const_iterator begin() const {return m_tank.begin();}
        const_iterator end() const {return m_tank.end();}

        SparseMatrix transpose() const;

        const Tank& tank() const {return m_tank;}

        Vector getlin(size_t i) const;

        void save(const char *filename) const;
        void load(const char *filename);

        void info() const;

        Vector operator*( const Vector &x ) const;
        Matrix operator*( const Matrix &m ) const;

    private:

        Tank m_tank;
    };

    inline Vector SparseMatrix::getlin(size_t i) const {
        assert(i<nlin());
        Vector v(ncol());
        for (size_t j=0;j<ncol();j++){
            const_iterator it = m_tank.find(std::make_pair<size_t, size_t>(i, j));
            if (it != m_tank.end()) v(j)=it->second;
            else v(j)=0.0;
        }
        return v;
    }
}
#endif  //! OPENMEEG_SPARSE_MATRIX_H
