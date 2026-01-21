// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <vector>

#include <vector.h>
#include <matrix.h>
#include <symmatrix.h>
#include <geometry.h>
#include <sensors.h>
#include <integrator.h>

#include <sparse_matrix.h>

/// @file
/// @brief Various helper functions for assembling matrices.

namespace OpenMEEG {

    // For ADAPT_LHS change the 0 in Integrator below into 10
    // It would be nice to define some constant integrators for the default values but swig does not like them.

    OPENMEEG_EXPORT SymMatrix HeadMat(const Geometry& geo,const Integrator& integrator=Integrator(3,0,0.005));
    OPENMEEG_EXPORT Matrix SurfSourceMat(const Geometry& geo,Mesh& sources,const Integrator& integrator=Integrator(3,0,0.005));

    OPENMEEG_EXPORT Matrix
    MonopoleSourceMat(const Geometry& geo,const Matrix& monopoles,const Integrator& integrator,const std::string& domain_name);
    OPENMEEG_EXPORT Matrix
    MonopoleSourceMat(const Geometry& geo,const Matrix& monopoles,const std::string& domain_name);

    OPENMEEG_EXPORT Matrix
    DipSourceMat(const Geometry& geo,const Matrix& dipoles,const Integrator& integrator,const std::string& domain_name);
    OPENMEEG_EXPORT Matrix
    DipSourceMat(const Geometry& geo,const Matrix& dipoles,const std::string& domain_name);

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
