/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
lastrevision      : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
*/

#include "matrice.h"

class HMEG_matrice : public virtual matrice
{
public:
    HMEG_matrice (const symmatrice& LhsInvMatrix,const matrice& RhsMatrix, const matrice& V2MegMatrix, const matrice& S2MegMatrix);
    virtual ~HMEG_matrice () {};
};

class HEEG_matrice : public virtual matrice
{
public:
    HEEG_matrice (const symmatrice& LhsInvMatrix,const matrice& RhsMatrix, const matrice& V2EegMatrix);
    virtual ~HEEG_matrice () {};
};

inline void assemble_gain_EEG(matrice& EegGainMatrix,const symmatrice& LhsInvMatrix,const matrice& RhsMatrix, const matrice& V2EegMatrix) {
    matrice reducedLhsInvMatrix = matrice(LhsInvMatrix)(0,LhsInvMatrix.nlin()-1,0,RhsMatrix.nlin()-1);
    EegGainMatrix = (V2EegMatrix*reducedLhsInvMatrix)*RhsMatrix;
}

inline void assemble_gain_MEG(matrice& MegGainMatrix,const symmatrice& LhsInvMatrix,const matrice& RhsMatrix, const matrice& V2MegMatrix, const matrice& S2MegMatrix) {
    matrice reducedLhsInvMatrix = matrice(LhsInvMatrix)(0,LhsInvMatrix.nlin()-1,0,RhsMatrix.nlin()-1);
    MegGainMatrix = S2MegMatrix+(V2MegMatrix*reducedLhsInvMatrix)*RhsMatrix;
}

HMEG_matrice::HMEG_matrice(const symmatrice& LhsInvMatrix,const matrice& RhsMatrix, const matrice& V2MegMatrix, const matrice& S2MegMatrix) {
    assemble_gain_MEG(*this,LhsInvMatrix,RhsMatrix,V2MegMatrix,S2MegMatrix);
}

HEEG_matrice::HEEG_matrice(const symmatrice& LhsInvMatrix,const matrice& RhsMatrix, const matrice& V2EegMatrix) {
    assemble_gain_EEG(*this,LhsInvMatrix,RhsMatrix,V2EegMatrix);
}
