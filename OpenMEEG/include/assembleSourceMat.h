// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include "matrix.h"
#include "geometry.h"
#include "integrator.h"
#include "sensors.h"

namespace OpenMEEG {

    // For ADAPT_LHS change the 0 in Integrator below into 10
    // It would be nice to define some constant integrators for the default values but swig does not like them.

    OPENMEEG_EXPORT Matrix SurfSourceMat(const Geometry& geo,Mesh& sources,const Integrator& integrator=Integrator(3,0,0.005));

    template <typename Source>
    OPENMEEG_EXPORT Matrix
    SourceMatrix(const Geometry& geo,const Matrix& sources,const Integrator& integrator,const std::string& domain_name);

    template <typename Source>
    Matrix SourceMatrix(const Geometry& geo,const Matrix& sources,const std::string& domain_name) {
        return SourceMatrix<Source>(geo,sources,Integrator(3,10,0.001),domain_name);
    }

    OPENMEEG_EXPORT Matrix EITSourceMat(const Geometry& geo,const Sensors& electrodes,const Integrator& integrator=Integrator(3,0,0.005));
    OPENMEEG_EXPORT Matrix DipSource2InternalPotMat(const Geometry& geo,const Matrix& dipoles,const Matrix& points,const std::string& domain_name="");
}
