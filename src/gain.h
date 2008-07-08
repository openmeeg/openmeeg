/* FILE: $Id$ */

/*
Project Name : OpenMEEG

version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

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

#include "matrix.h"
#include "sparse_matrix.h"

class HMEG_matrix : public virtual Matrix
{
public:
    HMEG_matrix (const SymMatrix& LhsInvMatrix,const Matrix& RhsMatrix, const Matrix& vToMEGMatrix, const Matrix& sToMEGMatrix);
    virtual ~HMEG_matrix () {};
};

class HEEG_matrix : public virtual Matrix
{
public:
    HEEG_matrix (const SymMatrix& LhsInvMatrix,const Matrix& RhsMatrix, const SparseMatrix& vToEEGMatrix);
    virtual ~HEEG_matrix () {};
};

inline void assemble_gain_EEG(Matrix& EEGGainMatrix,const SymMatrix& LhsInvMatrix,const Matrix& RhsMatrix, const SparseMatrix& vToEEGMatrix) {
    Matrix reducedLhsInvMatrix = LhsInvMatrix(0,LhsInvMatrix.nlin()-1,0,RhsMatrix.nlin()-1);
    EEGGainMatrix = (vToEEGMatrix*reducedLhsInvMatrix)*RhsMatrix;
}

inline void assemble_gain_MEG(Matrix& MEGGainMatrix,const SymMatrix& LhsInvMatrix,const Matrix& RhsMatrix, const Matrix& vToMEGMatrix, const Matrix& sToMEGMatrix) {
    Matrix reducedLhsInvMatrix = LhsInvMatrix(0,LhsInvMatrix.nlin()-1,0,RhsMatrix.nlin()-1);
    MEGGainMatrix = sToMEGMatrix+(vToMEGMatrix*reducedLhsInvMatrix)*RhsMatrix;
}

HMEG_matrix::HMEG_matrix(const SymMatrix& LhsInvMatrix,const Matrix& RhsMatrix, const Matrix& vToMEGMatrix, const Matrix& sToMEGMatrix) {
    assemble_gain_MEG(*this,LhsInvMatrix,RhsMatrix,vToMEGMatrix,sToMEGMatrix);
}

HEEG_matrix::HEEG_matrix(const SymMatrix& LhsInvMatrix,const Matrix& RhsMatrix, const SparseMatrix& vToEEGMatrix) {
    assemble_gain_EEG(*this,LhsInvMatrix,RhsMatrix,vToEEGMatrix);
}
