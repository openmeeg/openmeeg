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

#ifndef H_SPARSE_MATRICE
#define H_SPARSE_MATRICE

#include "sparse_matrice_dcl.h"

template <>
inline vecteur sparse_matrice::operator*(const vecteur &x) const
{
    vecteur ret(nlin());
    ret.set(0);

    for(idxType i=0;i<_nnzrows;i++)
    {
        MyCell<sparse_matrice::valType,sparse_matrice::idxType>* startRow=_RowEntry[_nzrows[i]];
        if(startRow!=0)
        {
            double total=0;
            while(startRow!=0)
            {
                total+=startRow->val*x(startRow->j);
                startRow=startRow->right;
            }
            ret(_nzrows[i])=total;
        }
    }

    return ret;
}

template <>
inline matrice sparse_matrice::operator*(const matrice &mat) const
{
    std::cout << "toto" << std::endl;
    std::cout << ncol() << std::endl;
    std::cout << mat.nlin() << std::endl;

    assert(ncol()==mat.nlin());
    matrice out(nlin(),mat.ncol());
    out.set(0.0);

    sparse_matrice_iterator it(*this);
    for(it.begin(); !it.end(); it.next()) {
        idxType ii = it.current()->i;
        idxType jj = it.current()->j;
        valType val = it.current()->val;
        for(size_t k = 0; k < mat.ncol(); ++k) {
            out(ii,k) += val * mat(jj,k);
        }
    }

    return out;
}

#endif
