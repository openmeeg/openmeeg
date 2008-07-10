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

class Head_matrix : public virtual SymMatrix
{
public:
    Head_matrix (const Geometry &geo, const int GaussOrder);
    virtual ~Head_matrix () {};
};

class SurfSource_matrix : public virtual Matrix
{
public:
    SurfSource_matrix (const Geometry &geo, const Mesh & sources, const int GaussOrder);
    virtual ~SurfSource_matrix () {};
};

class DipSource_matrix : public virtual Matrix
{
public:
    DipSource_matrix (const Geometry &geo, std::vector<Vect3> Rs, std::vector<Vect3> Qs, const int GaussOrder);
    virtual ~DipSource_matrix () {};
};

class DipSourceGrad_matrix : public virtual Matrix
{
public:
    DipSourceGrad_matrix (const Geometry &geo, std::vector<Vect3> Rs, std::vector<Vect3> Qs, const int GaussOrder);
    virtual ~DipSourceGrad_matrix () {};
};

class Surf2Vol_matrix : public virtual Matrix
{
public:
    Surf2Vol_matrix (const Geometry &geo, const Matrix &points);
    virtual ~Surf2Vol_matrix () {};
};
class Head2EEG_matrix : public virtual SparseMatrix
{
public:
    Head2EEG_matrix (const Geometry &geo, const Matrix& patches);
    virtual ~Head2EEG_matrix () {};
};

class Head2MEG_matrix : public virtual Matrix
{
public:
    Head2MEG_matrix (const Geometry &geo, const Sensors& sensors);
    virtual ~Head2MEG_matrix () {};
};

class SurfSource2MEG_matrix : public virtual Matrix
{
public:
    SurfSource2MEG_matrix (const Mesh& sources, const Sensors& sensors);
    virtual ~SurfSource2MEG_matrix () {};
};

class DipSource2MEG_matrix : public virtual Matrix
{
public:
    DipSource2MEG_matrix(const Matrix &dipoles, const Sensors &sensors);
    virtual ~DipSource2MEG_matrix () {};
};

void assemble_EITsource(const Geometry &geo, Matrix &mat, Matrix &airescalp, const int GaussOrder);

#endif /* ASSEMBLE_H */
