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
