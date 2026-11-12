// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <vector.h>
#include <matrix.h>
#include <operators.h>
#include <sensors.h>
#include <dipole.h>
#include <OMExceptions.H>
#include <constants.h>
#include <assembleSourceMat.h>

namespace OpenMEEG {

    Matrix SurfSourceMat(const Geometry& geo,Mesh& source_mesh,const Integrator& integrator) {

        // Check that there is no overlapping between the geometry and the source mesh.

        if (!geo.check(source_mesh))
            throw OverlappingSourceMesh();

        // By the previous test, source_mesh is included in a domain of the geometry.
        // Find this domain.

        const Domain& domain = geo.domain(*source_mesh.vertices().front());

        // Set it as an outermost (to tell operatorN it doesn't belong to the geometry).

        source_mesh.set_external_property(true);
        source_mesh.set_current_barrier_property(true);

        log_stream(INFORMATION) << std::endl
                                << "assemble SurfSourceMat with " << source_mesh.vertices().size()
                                << " source_mesh located in domain \"" << domain.name() << "\"." << std::endl
                                << std::endl;

        Matrix mat(geo.num_unknowns(),source_mesh.vertices().size());
        mat = 0.0;

        const double L  = -1.0/domain.conductivity();
        for (const auto& boundary : domain.boundaries()) {
            const double factorN = (boundary.is_outer_boundary()) ? K : -K;
            for (const auto& oriented_mesh : boundary.interface().oriented_meshes()) {
                const Mesh& mesh = oriented_mesh.mesh();

                NonDiagonalBlockWriter writer(mesh,source_mesh,integrator);

                // First block is nVertexFirstLayer x source_mesh.vertices().size()
                const double coeffN = factorN*oriented_mesh.orientation();
                writer.set_N_block(coeffN,mat);
                // Second block is nFacesFirstLayer x source_mesh.vertices().size()
                writer.D(coeffN*L,mat);
            }
        }

        return mat;
    }

    Matrix EITSourceMat(const Geometry& geo,const Sensors& electrodes,const Integrator& integrator) {

        // Matrix to be applied to the scalp-injected current to obtain the source term of the EIT foward problem,
        // following article BOUNDARY ELEMENT FORMULATION FOR ELECTRICAL IMPEDANCE TOMOGRAPHY
        // (eq.14 (do not look at eq.16 since there is a mistake: D_23 -> S_23))
        // rhs = [0 ... 0  -D*_23  sigma_3^(-1)S_23  -I_33/2.+D*_33]

        SymMatrix transmat(geo.num_parameters()); // TODO: num_unknowns instead ?
        transmat = 0.0;

        // This is an overkill. Can we limit the computation only to injection triangles ?
        // We use only the few lines that correspond to injection triangles.

        for (const Mesh* mp : geo.active_meshes()) {
            const Mesh& mesh = *mp;
            if (mesh.is_current_barrier()) {
                const DiagonalBlockWriter writer(mesh,integrator);
                writer.set_D_block(transmat);      // D*_33 part of the formula.
                writer.addIdentity(-0.5,transmat); // I_33 part of the formula.
            }
        }

        for (const auto& mp : geo.interacting_meshes()) {
            const Mesh& mesh1 = mp(0);
            const Mesh& mesh2 = mp(1);

            if (mesh1.is_current_barrier()) {
                const ConductivityCoefficients cond_coeffs(mp);
                const NonDiagonalBlockWriter writer(mesh1,mesh2,integrator);
                writer.D(-cond_coeffs.D(),transmat); // D*_23 part of the formula.
                writer.S(-cond_coeffs.S(),transmat); // S_23 part of the formula.
            }
        }

        const size_t n_sensors = electrodes.getNumberOfSensors();
        Matrix mat(geo.num_unknowns(),n_sensors);
        mat = 0.0;

        for (size_t ielec=0; ielec<n_sensors; ++ielec)
            for (const auto& triangle : electrodes.getInjectionTriangles(ielec)) {
                // To ensure exactly no accumulation of currents. w = electrode_area/triangle_area (~= 1)
                // If no radius is given, we assume the user wants to specify an intensity not a density of current.
                const double coeff = (almost_equal(electrodes.getRadii()(ielec),0.0)) ? 1.0/triangle.area() : electrodes.getWeights()(ielec);
                for (size_t i=0; i<mat.nlin(); ++i)
                    mat(i,ielec) += transmat(triangle.index(),i)*coeff;
            }
        return mat;
    }

    Matrix DipSource2InternalPotMat(const Geometry& geo,const Matrix& dipoles,const Matrix& points,const std::string& domain_name) {

        // Points with one more column for the index of the domain they belong

        std::vector<const Domain*> points_domain;
        std::vector<Vect3>         pts;
        for (unsigned i=0; i<points.nlin(); ++i) {
            const Vect3   point(points(i,0),points(i,1),points(i,2));
            const Domain& domain = geo.domain(point);
            if (domain.conductivity()!=0.0) {
                points_domain.push_back(&domain);
                pts.push_back(point);
            } else {
                std::cerr << " DipSource2InternalPot: Point [ " << points.row(i)
                          << "] is outside the head. Point is dropped." << std::endl;
            }
        }

        Matrix mat(pts.size(),dipoles.nlin());
        mat = 0.0;

        for (unsigned iDIP=0; iDIP<dipoles.nlin(); ++iDIP) {
            const Dipole dipole(iDIP,dipoles);

            const Domain& domain = (domain_name=="") ? geo.domain(dipole.position()) : geo.domain(domain_name);
            const double  coeff  = K/domain.conductivity();

            for (unsigned iPTS=0; iPTS<pts.size(); ++iPTS)
                if (points_domain[iPTS]==&domain)
                    mat(iPTS,iDIP) += coeff*dipole.potential(pts[iPTS]);
        }
        return mat;
    }
}
