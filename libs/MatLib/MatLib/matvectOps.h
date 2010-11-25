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

#ifndef OPENMEEG_MATVECTOPS_H
#define OPENMEEG_MATVECTOPS_H

#include "vector.h"
#include "matrix.h"
#include "symmatrix.h"
#include "sparse_matrix.h"
#include "diagmatrix.h"

namespace OpenMEEG {

    inline SymMatrix SymMatrix::operator*(const SymMatrix &m) const {
        assert(nlin()==m.nlin());
    #ifdef HAVE_BLAS
        Matrix D(*this);
        Matrix B(m);
        Matrix C(nlin(),nlin());
        DSYMM(CblasLeft,CblasUpper,(int)nlin(),(int)B.ncol(),1.,D.data(),(int)D.ncol(),B.data(),(int)B.nlin(),0,C.data(),(int)C.nlin());
        return SymMatrix(C);
    #else
        SymMatrix C(nlin());
        for (size_t j=0;j<m.ncol();j++){
            for (size_t i=0;i<ncol();i++)
            {
                C(i,j)=0;
                for (size_t k=0;k<ncol();k++)
                    C(i,j)+=(*this)(i,k)*m(k,j);
            }
        }
        return C;
    #endif
    }

    inline Matrix SymMatrix::operator*(const Matrix &B) const {
        assert(ncol()==B.nlin());
        Matrix C(nlin(),B.ncol());
    #ifdef HAVE_BLAS
        Matrix D(*this);
        DSYMM(CblasLeft,CblasUpper ,(int)nlin(), (int)B.ncol(), 1. , D.data(), (int)D.ncol(), B.data(), (int)B.nlin(), 0, C.data(),(int)C.nlin());
    #else
        for (size_t j=0;j<B.ncol();j++){
            for (size_t i=0;i<ncol();i++)
            {
                C(i,j)=0;
                for (size_t k=0;k<ncol();k++)
                    C(i,j)+=(*this)(i,k)*B(k,j);
            }
        }
    #endif
        return C;
    }

    inline Matrix SymMatrix::solve(Matrix &RHS) const {
    #ifdef HAVE_LAPACK
        SymMatrix A(*this,DEEP_COPY);
        // LU
        int *pivots=new int[nlin()];
        int Info;
        DSPTRF('U',A.nlin(),A.data(),pivots,Info);
        // Solve the linear system AX=B
        DSPTRS('U',A.nlin(),RHS.ncol(),A.data(),pivots,RHS.data(),A.nlin(),Info);
        return RHS;
    #else
        std::cerr << "!!!!! Not implemented : Try a GMres !!!!!" << std::endl;
        exit(1);
    #endif
    }
}
#endif // ! OPENMEEG_MATVECTOPS_H
