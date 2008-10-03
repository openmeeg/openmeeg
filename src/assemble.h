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

class OPENMEEG_EXPORT HeadMat : public virtual SymMatrix
{
public:
    HeadMat (const Geometry &geo, const int GaussOrder);
    virtual ~HeadMat () {};
};

class OPENMEEG_EXPORT SurfSourceMat : public virtual Matrix
{
public:
    SurfSourceMat (const Geometry &geo, const Mesh & sources, const int GaussOrder);
    virtual ~SurfSourceMat () {};
};

class OPENMEEG_EXPORT DipSourceMat : public virtual Matrix
{
public:
    DipSourceMat (const Geometry &geo, const Matrix& dipoles, const int GaussOrder);
    virtual ~DipSourceMat () {};
};

class OPENMEEG_EXPORT DipSourceGradMat : public virtual Matrix
{
public:
    DipSourceGradMat (const Geometry &geo, std::vector<Vect3> Rs, std::vector<Vect3> Qs, const int GaussOrder);
    virtual ~DipSourceGradMat () {};
};

class OPENMEEG_EXPORT Surf2VolMat : public virtual Matrix
{
public:
    Surf2VolMat (const Geometry &geo, const Matrix &points);
    virtual ~Surf2VolMat () {};
};

class OPENMEEG_EXPORT Head2EEGMat : public virtual SparseMatrix
{
public:
    Head2EEGMat (const Geometry &geo, const Matrix& patches);
    virtual ~Head2EEGMat () {};
};

class OPENMEEG_EXPORT Head2MEGMat : public virtual Matrix
{
public:
    Head2MEGMat (const Geometry &geo, const Sensors& sensors);
    virtual ~Head2MEGMat () {};
};

class OPENMEEG_EXPORT SurfSource2MEGMat : public virtual Matrix
{
public:
    SurfSource2MEGMat (const Mesh& sources, const Sensors& sensors);
    virtual ~SurfSource2MEGMat () {};
};

class OPENMEEG_EXPORT DipSource2MEGMat : public virtual Matrix
{
public:
    DipSource2MEGMat(const Matrix &dipoles, const Sensors &sensors);
    virtual ~DipSource2MEGMat () {};
};

OPENMEEG_EXPORT void assemble_EITsource(const Geometry &geo, Matrix &mat, Matrix &airescalp, const int GaussOrder);

#endif /* ASSEMBLE_H */
