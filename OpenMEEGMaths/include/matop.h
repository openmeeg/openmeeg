// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include "matrix.h"
#include "sparse_matrix.h"

namespace OpenMEEG {

    inline Matrix
    nullspace_projector(const Matrix& M) {
        const size_t Nl = M.nlin();
        const size_t Nc = M.ncol();

        Matrix U,V;
        SparseMatrix D;
        M.svd(U,D,V);

        // Set S to 0 everywhere, except in the last part of the diag:

        SparseMatrix S(Nc,Nc);
        for (unsigned i=Nl;i<Nc;++i)
            S(i,i) = 1.0;

        return (V.transpose()*S)*V; // P is a projector: P^2 = P and mat*P*X = 0
    }
}
