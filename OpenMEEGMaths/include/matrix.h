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

#include <OpenMEEGMathsConfig.h>

#include <linop.h>

namespace OpenMEEG {

    class OPENMEEGMATHS_EXPORT Matrix: public LinOp {
    protected:

        // friend class Vector;

        // utils::RCPtr<LinOpValue> value;

        // explicit Matrix(const Matrix& A,const size_t M): LinOp(A.nlin(),M,FULL,2),value(A.value) { }

    public:
        Matrix(): LinOp(0,0,FULL,2) { }

        // Matrix(): LinOp(0,0,FULL,2),value() { }
        // Matrix(const char* fname): LinOp(0,0,FULL,2),value() { this->load(fname); }
        // Matrix(const size_t M,const size_t N): LinOp(M,N,FULL,2),value(new LinOpValue(N*M)) { }
        // Matrix(const Matrix& A,const DeepCopy): LinOp(A.nlin(),A.ncol(),FULL,2),value(new LinOpValue(A.size(),A.data())) { }

        // bool empty() const { return value->empty(); }

        // size_t size() const { return nlin()*ncol(); };

        // double* data() const { return value->data; }
        double* data() const { return 0; }

        // const Matrix& set(const double d);
        // void load(const char *filename);
        // void save(const char *filename) const;

    };

}
