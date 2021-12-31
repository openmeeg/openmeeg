/*
Project Name: OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

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

#pragma once

#include <vector>

#include <vector.h>
#include <matrix.h>
#include <symmatrix.h>
#include <geometry.h>
#include <sensors.h>
#include <integrator.h>

#include <sparse_matrix.h>

namespace OpenMEEG {

    // For ADAPT_LHS change the 0 in Integrator below into 10
    // It would be nice to define some constant integrators for the default values but swig does not like them. 

    OPENMEEG_EXPORT SymMatrix HeadMat(const Geometry& geo,const Integrator& integrator=Integrator(3,0,0.005));
    OPENMEEG_EXPORT Matrix SurfSourceMat(const Geometry& geo,Mesh& sources,const Integrator& integrator=Integrator(3,0,0.005));

    OPENMEEG_EXPORT Matrix
    DipSourceMat(const Geometry& geo,const Matrix& dipoles,const Integrator& integrator=Integrator(3,10,0.001),const std::string& domain_name="");

    OPENMEEG_EXPORT Matrix EITSourceMat(const Geometry& geo,const Sensors& electrodes,const Integrator& integrator=Integrator(3,0,0.005));

    OPENMEEG_EXPORT Matrix Surf2VolMat(const Geometry& geo,const Matrix& points);

    OPENMEEG_EXPORT SparseMatrix Head2EEGMat(const Geometry& geo,const Sensors& electrodes);
    OPENMEEG_EXPORT SparseMatrix Head2ECoGMat(const Geometry& geo,const Sensors& electrodes,const Interface& i);

    inline SparseMatrix
    Head2ECoGMat(const Geometry& geo,const Sensors& electrodes,const std::string& id) { // Mainly for SWIG
        return Head2ECoGMat(geo,electrodes,geo.interface(id));
    }

    OPENMEEG_EXPORT Matrix Head2MEGMat(const Geometry& geo,const Sensors& sensors);
    OPENMEEG_EXPORT Matrix SurfSource2MEGMat(const Mesh& sources,const Sensors& sensors);
    OPENMEEG_EXPORT Matrix DipSource2MEGMat(const Matrix& dipoles,const Sensors& sensors);
    OPENMEEG_EXPORT Matrix DipSource2InternalPotMat(const Geometry& geo,const Matrix& dipoles,const Matrix& points,const std::string& domain_name="");

    OPENMEEG_EXPORT Matrix CorticalMat(const Geometry& geo,const SparseMatrix& M,const std::string& domain_name="CORTEX",
                                       const double alpha=-1.0,const double beta=-1.0,const std::string &filename="",
                                       const Integrator& integrator=Integrator(3,0,0.005));

    OPENMEEG_EXPORT Matrix CorticalMat2(const Geometry& geo,const SparseMatrix& M,const std::string& domain_name="CORTEX",
                                        const double gamma=1.0,const std::string &filename="",const Integrator& integrator=Integrator(3,0,0.005));
}
