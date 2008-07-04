/* FILE: $Id$ */

/*
Project Name : OpenMEEG

author            : $Author$
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

#include "matrice.h"
#include "sparse_matrice.h"

class HMEG_matrice : public virtual matrice
{
public:
    HMEG_matrice (const symmatrice& LhsInvMatrix,const matrice& RhsMatrix, const matrice& vToMEGMatrix, const matrice& sToMEGMatrix);
    virtual ~HMEG_matrice () {};
};

class HEEG_matrice : public virtual matrice
{
public:
    HEEG_matrice (const symmatrice& LhsInvMatrix,const matrice& RhsMatrix, const sparse_matrice& vToEEGMatrix);
    virtual ~HEEG_matrice () {};
};

inline void assemble_gain_EEG(matrice& EEGGainMatrix,const symmatrice& LhsInvMatrix,const matrice& RhsMatrix, const sparse_matrice& vToEEGMatrix) {
    matrice reducedLhsInvMatrix = LhsInvMatrix(0,LhsInvMatrix.nlin()-1,0,RhsMatrix.nlin()-1);
    EEGGainMatrix = (vToEEGMatrix*reducedLhsInvMatrix)*RhsMatrix;
}

inline void assemble_gain_MEG(matrice& MEGGainMatrix,const symmatrice& LhsInvMatrix,const matrice& RhsMatrix, const matrice& vToMEGMatrix, const matrice& sToMEGMatrix) {
    matrice reducedLhsInvMatrix = LhsInvMatrix(0,LhsInvMatrix.nlin()-1,0,RhsMatrix.nlin()-1);
    MEGGainMatrix = sToMEGMatrix+(vToMEGMatrix*reducedLhsInvMatrix)*RhsMatrix;
}

HMEG_matrice::HMEG_matrice(const symmatrice& LhsInvMatrix,const matrice& RhsMatrix, const matrice& vToMEGMatrix, const matrice& sToMEGMatrix) {
    assemble_gain_MEG(*this,LhsInvMatrix,RhsMatrix,vToMEGMatrix,sToMEGMatrix);
}

HEEG_matrice::HEEG_matrice(const symmatrice& LhsInvMatrix,const matrice& RhsMatrix, const sparse_matrice& vToEEGMatrix) {
    assemble_gain_EEG(*this,LhsInvMatrix,RhsMatrix,vToEEGMatrix);
}
