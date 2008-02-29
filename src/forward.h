/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
*/

#include "matrice.h"
#include "symmatrice.h"
#include "vecteur.h"

class Forward_matrice : public virtual matrice
{
public:
    Forward_matrice (const matrice& GainMatrix, const matrice& RealSourcesData, double NoiseLevel);
    virtual ~Forward_matrice () {};
};

void compute_forward(matrice& SimulatedData, const matrice& GainMatrix, const matrice& RealSourcesData, double NoiseLevel) {

    SimulatedData = GainMatrix * RealSourcesData;

    int nT = RealSourcesData.ncol();
    for(int frame=0;frame<nT;frame++)
    {
        for(size_t i=0;i<SimulatedData.nlin();i++) SimulatedData(i,frame) += NoiseLevel * gaussian();
    }
}

Forward_matrice::Forward_matrice(const matrice& GainMatrix, const matrice& RealSourcesData, double NoiseLevel) {
    compute_forward(*this,GainMatrix,RealSourcesData,NoiseLevel);
}