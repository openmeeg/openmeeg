// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <vector.h>
#include <matrix.h>
#include <danielsson.h>
#include <operators.h>
#include <assemble.h>
#include <sensors.h>
#include <OMExceptions.H>

#include <constants.h>

namespace OpenMEEG {

    Matrix SurfSourceMat(const Geometry& geo,Mesh& source_mesh,const Integrator& integrator) {

        // Check that there is no overlapping between the geometry and the source mesh.

        if (!geo.check(source_mesh))
            throw OverlappingSourceMesh();

        // The mesh is included in a domain of the geometry.

        const Domain& domain = geo.domain(*source_mesh.vertices().front());

        // Set it as an outermost (to tell _operarorN it doesn't belong to the geometry).

        source_mesh.outermost()       = true;
        source_mesh.current_barrier() = true;

        log_stream(INFORMATION) << std::endl
                                << "assemble SurfSourceMat with " << source_mesh.vertices().size()
                                << " source_mesh located in domain \"" << domain.name() << "\"." << std::endl
                                << std::endl;

        Matrix mat(geo.nb_parameters()-geo.nb_current_barrier_triangles(),source_mesh.vertices().size());
        mat.set(0.0);

        const double L  = -1.0/domain.conductivity();
        for (const auto& boundary : domain.boundaries()) {
            const double factorN = (boundary.inside()) ? K : -K;
            for (const auto& oriented_mesh : boundary.interface().oriented_meshes()) {
                const Mesh& mesh = oriented_mesh.mesh();

                NonDiagonalBlock operators(mesh,source_mesh,integrator);

                // First block is nVertexFistLayer*source_mesh.vertices().size()
                const double coeffN = factorN*oriented_mesh.orientation();
                operators.set_N_block(coeffN,mat);
                // Second block is nFacesFistLayer*source_mesh.vertices().size()
                operators.D(coeffN*L,mat);
            }
        }

        return mat;
    }

    Matrix
    MonopoleSourceMat(const Geometry& geo,const Matrix& monopoles,const Integrator& integrator,const std::string& domain_name) {

        const size_t size      = geo.nb_parameters()-geo.nb_current_barrier_triangles();
        const size_t n_monopoles = monopoles.nlin();

        Matrix rhs(size,n_monopoles);
        rhs.set(0.0);

        ProgressBar pb(n_monopoles);
        Vector rhs_col(rhs.nlin());
        for (unsigned s=0; s<n_monopoles; ++s,++pb) {
            const Monopole monopole(s,monopoles);
            const Domain domain = (domain_name=="") ? geo.domain(monopole.position()) : geo.domain(domain_name);

            //  Only consider monopoles in non-zero conductivity domain.

            const double cond = domain.conductivity();
            if (cond!=0.0) {
                rhs_col.set(0.0);
                for (const auto& boundary : domain.boundaries()) {
                    const double factorD = (boundary.inside()) ? K : -K;
                    for (const auto& oriented_mesh : boundary.interface().oriented_meshes()) {
                        //  Process the mesh.
                        const double coeffD = factorD*oriented_mesh.orientation();
                        const Mesh&  mesh   = oriented_mesh.mesh();
                        operatorMonopolePotDer(monopole,mesh,rhs_col,coeffD,integrator);

                        if (!oriented_mesh.mesh().current_barrier()) {
                            const double coeff = -coeffD/cond;;
                            operatorMonopolePot(monopole,mesh,rhs_col,coeff,integrator);
                        }
                    }
                }
                rhs.setcol(s,rhs_col);
            }
        }
        return rhs;
    }

    Matrix
    MonopoleSourceMat(const Geometry& geo,const Matrix& monopoles,const std::string& domain_name) {
        return MonopoleSourceMat(geo,monopoles,Integrator(3,10,0.001),domain_name);
    }

    Matrix
    DipSourceMat(const Geometry& geo,const Matrix& dipoles,const Integrator& integrator,const std::string& domain_name) {

        const size_t size      = geo.nb_parameters()-geo.nb_current_barrier_triangles();
        const size_t n_dipoles = dipoles.nlin();

        Matrix rhs(size,n_dipoles);
        rhs.set(0.0);

        ProgressBar pb(n_dipoles);
        Vector rhs_col(rhs.nlin());
        for (unsigned s=0; s<n_dipoles; ++s,++pb) {
            const Dipole dipole(s,dipoles);
            const Domain domain = (domain_name=="") ? geo.domain(dipole.position()) : geo.domain(domain_name);

            //  Only consider dipoles in non-zero conductivity domain.

            const double cond = domain.conductivity();
            if (cond!=0.0) {
                rhs_col.set(0.0);
                for (const auto& boundary : domain.boundaries()) {
                    const double factorD = (boundary.inside()) ? K : -K;
                    for (const auto& oriented_mesh : boundary.interface().oriented_meshes()) {
                        //  Process the mesh.
                        const double coeffD = factorD*oriented_mesh.orientation();
                        const Mesh&  mesh   = oriented_mesh.mesh();
                        operatorDipolePotDer(dipole,mesh,rhs_col,coeffD,integrator);

                        if (!oriented_mesh.mesh().current_barrier()) {
                            const double coeff = -coeffD/cond;;
                            operatorDipolePot(dipole,mesh,rhs_col,coeff,integrator);
                        }
                    }
                }
                rhs.setcol(s,rhs_col);
            }
        }
        return rhs;
    }

    Matrix
    DipSourceMat(const Geometry& geo,const Matrix& dipoles,const std::string& domain_name) {
        return DipSourceMat(geo,dipoles,Integrator(3,10,0.001),domain_name);
    }

    Matrix EITSourceMat(const Geometry& geo,const Sensors& electrodes,const Integrator& integrator) {

        // Matrix to be applied to the scalp-injected current to obtain the source term of the EIT foward problem,
        // following article BOUNDARY ELEMENT FORMULATION FOR ELECTRICAL IMPEDANCE TOMOGRAPHY
        // (eq.14 (do not look at eq.16 since there is a mistake: D_23 -> S_23))
        // rhs = [0 ... 0  -D*_23  sigma_3^(-1)S_23  -I_33/2.+D*_33]

        SymMatrix transmat(geo.nb_parameters());
        transmat.set(0.0);

        // This is an overkill. Can we limit the computation only to injection triangles ?
        // We use only the few lines that correspond to injection triangles.

        for (const auto& mp : geo.communicating_mesh_pairs()) {
            const Mesh& mesh1 = mp(0);
            const Mesh& mesh2 = mp(1);

            if (mesh1.current_barrier()) {
                const NonDiagonalBlock operators(mesh1,mesh2,integrator);
                const int orientation = geo.relative_orientation(mesh1,mesh2);
                operators.D(K*orientation,transmat); // D23 or D33 of the formula.
                if (&mesh1==&mesh2) { // I_33 of the formula, orientation is necessarily 1.
                    DiagonalBlock block(mesh1,integrator);
                    block.addIdentity(-0.5,transmat);
                } else { // S_2 of the formula.
                    operators.S(-K*orientation*geo.sigma_inv(mesh1,mesh2),transmat);
                }
            }
        }

        const size_t n_sensors = electrodes.getNumberOfSensors();
        Matrix mat(geo.nb_parameters()-geo.nb_current_barrier_triangles(),n_sensors);
        mat.set(0.0);

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
                std::cerr << " DipSource2InternalPot: Point [ " << points.getlin(i)
                          << "] is outside the head. Point is dropped." << std::endl;
            }
        }

        Matrix mat(pts.size(),dipoles.nlin());
        mat.set(0.0);

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
