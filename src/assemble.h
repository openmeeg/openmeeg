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

#ifndef _ASSEMBLE_H_
#define _ASSEMBLE_H_

#include <vector>
using namespace std;

#include "vecteur.h"
#include "matrice.h"
#include "symmatrice.h"
#include "geometry.h"
#include "sensors.h"

class LHS_matrice : public virtual symmatrice
{
public:
    LHS_matrice (const Geometry &geo, const int GaussOrder);
    virtual ~LHS_matrice () {};
};

class RHS_matrice : public virtual matrice
{
public:
    RHS_matrice (const Geometry &geo, const Mesh & sources, const int GaussOrder);
    virtual ~RHS_matrice () {};
};

class RHSdip_matrice : public virtual matrice
{
public:
    RHSdip_matrice (const Geometry &geo, vector<Vect3> Rs, vector<Vect3> Qs, const int GaussOrder);
    virtual ~RHSdip_matrice () {};
};

class RHSdip_grad_matrice : public virtual matrice
{
public:
    RHSdip_grad_matrice (const Geometry &geo, vector<Vect3> Rs, vector<Vect3> Qs, const int GaussOrder);
    virtual ~RHSdip_grad_matrice () {};
};

class SurfToVol_matrice : public virtual matrice
{
public:
    SurfToVol_matrice (const Geometry &geo, const matrice &points);
    virtual ~SurfToVol_matrice () {};
};
class vToEEG_matrice : public virtual matrice
{
public:
    vToEEG_matrice (const Geometry &geo, const matrice& patches);
    virtual ~vToEEG_matrice () {};
};

class vToMEG_matrice : public virtual matrice
{
public:
    vToMEG_matrice (const Geometry &geo, const Sensors& sensors);
    virtual ~vToMEG_matrice () {};
};

class sToMEG_matrice : public virtual matrice
{
public:
    sToMEG_matrice (const Mesh& sources, const Sensors& sensors);
    virtual ~sToMEG_matrice () {};
};

class sToMEGdip_matrice : public virtual matrice
{
public:
    sToMEGdip_matrice(const matrice &dipoles, const Sensors &sensors);
    virtual ~sToMEGdip_matrice () {};
};

void assemble_EITsource(const Geometry &geo, matrice &mat, matrice &airescalp, const int GaussOrder);
#endif /* _ASSEMBLE_H_ */
