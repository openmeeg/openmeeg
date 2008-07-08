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

#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include <vector>

#include "vector.h"
#include "matrix.h"
#include "symmatrix.h"
#include "geometry.h"
#include "sensors.h"

class LHS_matrix : public virtual SymMatrix
{
public:
    LHS_matrix (const Geometry &geo, const int GaussOrder);
    virtual ~LHS_matrix () {};
};

class RHS_matrix : public virtual Matrix
{
public:
    RHS_matrix (const Geometry &geo, const Mesh & sources, const int GaussOrder);
    virtual ~RHS_matrix () {};
};

class RHSdip_matrix : public virtual Matrix
{
public:
    RHSdip_matrix (const Geometry &geo, std::vector<Vect3> Rs, std::vector<Vect3> Qs, const int GaussOrder);
    virtual ~RHSdip_matrix () {};
};

class RHSdip_grad_matrix : public virtual Matrix
{
public:
    RHSdip_grad_matrix (const Geometry &geo, std::vector<Vect3> Rs, std::vector<Vect3> Qs, const int GaussOrder);
    virtual ~RHSdip_grad_matrix () {};
};

class SurfToVol_matrix : public virtual Matrix
{
public:
    SurfToVol_matrix (const Geometry &geo, const Matrix &points);
    virtual ~SurfToVol_matrix () {};
};
class vToEEG_matrix : public virtual SparseMatrix
{
public:
    vToEEG_matrix (const Geometry &geo, const Matrix& patches);
    virtual ~vToEEG_matrix () {};
};

class vToMEG_matrix : public virtual Matrix
{
public:
    vToMEG_matrix (const Geometry &geo, const Sensors& sensors);
    virtual ~vToMEG_matrix () {};
};

class sToMEG_matrix : public virtual Matrix
{
public:
    sToMEG_matrix (const Mesh& sources, const Sensors& sensors);
    virtual ~sToMEG_matrix () {};
};

class sToMEGdip_matrix : public virtual Matrix
{
public:
    sToMEGdip_matrix(const Matrix &dipoles, const Sensors &sensors);
    virtual ~sToMEGdip_matrix () {};
};

void assemble_EITsource(const Geometry &geo, Matrix &mat, Matrix &airescalp, const int GaussOrder);

#endif /* ASSEMBLE_H */
