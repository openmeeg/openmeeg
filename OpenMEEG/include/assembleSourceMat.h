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
#include "operators.h"
#include "progressbar.h"

#include <constants.h>

namespace OpenMEEG {

    // For ADAPT_LHS change the 0 in Integrator below into 10
    // It would be nice to define some constant integrators for the default values but swig does not like them.

    OPENMEEG_EXPORT Matrix SurfSourceMat(const Geometry& geo,Mesh& sources,const Integrator& integrator=Integrator(3,0,0.005));

    template <typename Source>
    Matrix
    SourceMatrix(const Geometry& geo,const Matrix& sources,const Integrator& integrator,const std::string& domain_name) {

        const size_t size      = geo.nb_parameters()-geo.nb_current_barrier_triangles();
        const size_t n_sources = sources.nlin();

        Matrix rhs(size,n_sources);
        rhs = 0.0;

        ProgressBar pb(n_sources);
        for (unsigned s=0; s<n_sources; ++s,++pb) {
            const Source source(s,sources);
            const Domain domain = (domain_name=="") ? geo.domain(source.position()) : geo.domain(domain_name);

            //  Only consider sources in non-zero conductivity domains.

            const double cond = domain.conductivity();
            if (cond!=0.0) {
                for (const auto& boundary : domain.boundaries()) {
                    const double factorD = (boundary.inside()) ? -K : K;
                    for (const auto& oriented_mesh : boundary.interface().oriented_meshes()) {
                        //  Process the mesh.
                        if (s==4)
                            std::cerr << "Before: " << rhs.column(s) << std::endl;
                        const double coeffD = factorD*oriented_mesh.orientation();
                        const Mesh&  mesh   = oriented_mesh.mesh();
                        rhs.column(s) += coeffD*operatorPotentialDerivative(size,source,mesh,integrator);

                        if (s==4)
                            std::cerr << "Middle: " << rhs.column(s) << std::endl;

                        if (!oriented_mesh.mesh().current_barrier()) {
                            const double coeff = coeffD/cond;;
                            rhs.column(s) += coeff*operatorPotential(size,source,mesh,integrator);
                        }

                        if (s==4)
                            std::cerr << "After: " << rhs.column(s) << std::endl;
                        if (s==4)
                            std::cerr << "=======================================================" << std::endl;
                    }
                }
            }
        }
        return rhs;
    }

    template <typename Source>
    Matrix SourceMatrix(const Geometry& geo,const Matrix& sources,const std::string& domain_name) {
        return SourceMatrix<Source>(geo,sources,Integrator(3,10,0.001),domain_name);
    }

    OPENMEEG_EXPORT Matrix EITSourceMat(const Geometry& geo,const Sensors& electrodes,const Integrator& integrator=Integrator(3,0,0.005));
    OPENMEEG_EXPORT Matrix DipSource2InternalPotMat(const Geometry& geo,const Matrix& dipoles,const Matrix& points,const std::string& domain_name="");
}
