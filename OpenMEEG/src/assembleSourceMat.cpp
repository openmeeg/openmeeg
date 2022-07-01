/*
Project Name : OpenMEEG

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

#include <vector.h>
#include <matrix.h>
#include <danielsson.h>
#include <operators.h>
#include <assemble.h>
#include <sensors.h>

#include <constants.h>

namespace OpenMEEG {

    SurfSourceMat::SurfSourceMat(const Geometry& geo,Mesh& source_mesh,const unsigned gauss_order) {

        Matrix& mat = *this;

        mat = Matrix((geo.nb_parameters()-geo.nb_current_barrier_triangles()),source_mesh.vertices().size());
        mat.set(0.0);

        // Check that there is no overlapping between the geometry and the source mesh.

        if (!geo.check(source_mesh)) {
            std::cerr << "Error: source mesh overlapps the geometry" << std::endl;
            return;
        }

        // The mesh is included in a domain of the geometry.

        const Domain& domain = geo.domain(*source_mesh.vertices().front());

        // Set it as an outermost (to tell _operarorN it doesn't belong to the geometry).

        source_mesh.outermost()       = true;
        source_mesh.current_barrier() = true;

        std::cout << std::endl
                  << "assemble SurfSourceMat with " << source_mesh.vertices().size()
                  << " source_mesh located in domain \"" << domain.name() << "\"." << std::endl
                  << std::endl;

        const double K  = 1.0/(4*Pi);
        const double L  = -1.0/domain.conductivity();
        for (const auto& boundary : domain.boundaries()) {
            const double factorN = (boundary.inside()) ? K : -K;
            for (const auto& oriented_mesh : boundary.interface().oriented_meshes()) {
                const Mesh& mesh = oriented_mesh.mesh();
                // First block is nVertexFistLayer*source_mesh.vertices().size()
                const double coeffN = factorN*oriented_mesh.orientation();
                operatorN(mesh,source_mesh,mat,coeffN,gauss_order);
                // Second block is nFacesFistLayer*source_mesh.vertices().size()
                operatorD(mesh,source_mesh,mat,coeffN*L,gauss_order);
            }
        }
    }

    DipSourceMat::DipSourceMat(const Geometry& geo,const Matrix& dipoles,const unsigned gauss_order,
                               const bool adapt_rhs,const std::string& domain_name)
    {
        Matrix& rhs = *this;

        const size_t size      = geo.nb_parameters()-geo.nb_current_barrier_triangles();
        const size_t n_dipoles = dipoles.nlin();

        rhs = Matrix(size,n_dipoles);
        rhs.set(0.0);

        ProgressBar pb(n_dipoles);
        Vector rhs_col(rhs.nlin());
        for (unsigned s=0; s<n_dipoles; ++s,++pb) {
            const Vect3 r(dipoles(s,0),dipoles(s,1),dipoles(s,2));
            const Vect3 q(dipoles(s,3),dipoles(s,4),dipoles(s,5));

            const Domain domain = (domain_name=="") ? geo.domain(r) : geo.domain(domain_name);

            //  Only consider dipoles in non-zero conductivity domain.

            const double cond = domain.conductivity();
            if (cond!=0.0) {
                rhs_col.set(0.0);
                const double K = 1.0/(4*Pi);
                for (const auto& boundary : domain.boundaries()) { //  Iterate over the domain's interfaces (half-spaces)
                    const double factorD = (boundary.inside()) ? K : -K;
                    for (const auto& oriented_mesh : boundary.interface().oriented_meshes()) { //  Iterate over the meshes of the interface
                        //  Treat the mesh.
                        const double coeffD = factorD*oriented_mesh.orientation();
                        const Mesh&  mesh   = oriented_mesh.mesh();
                        operatorDipolePotDer(r,q,mesh,rhs_col,coeffD,gauss_order,adapt_rhs);

                        if (!oriented_mesh.mesh().current_barrier()) {
                            const double coeff = -coeffD/cond;;
                            operatorDipolePot(r,q,mesh,rhs_col,coeff,gauss_order,adapt_rhs);
                        }
                    }
                }
                rhs.setcol(s,rhs_col);
            }
        }
    }

    EITSourceMat::EITSourceMat(const Geometry& geo,const Sensors& electrodes,const unsigned gauss_order) {
        Matrix& mat = *this;

        //  A Matrix to be applied to the scalp-injected current to obtain the Source Term of the EIT foward problem.
        // following article BOUNDARY ELEMENT FORMULATION FOR ELECTRICAL IMPEDANCE TOMOGRAPHY
        // (eq.14 (do not look at eq.16 since there is a mistake: D_23 -> S_23))
        // rhs = [0 ... 0  -D*_23  sigma_3^(-1)S_23  -I_33/2.+D*_33]

        size_t n_sensors = electrodes.getNumberOfSensors();

        const double K = 1.0/(4*Pi);

        //  transmat = a big SymMatrix of which mat = part of its transpose.
        SymMatrix transmat(geo.nb_parameters());
        transmat.set(0.0);
        mat = Matrix((geo.nb_parameters()-geo.nb_current_barrier_triangles()), n_sensors);
        mat.set(0.0);

        for (const auto& mesh1 : geo.meshes())
            if (mesh1.current_barrier())
                for (const auto& mesh2 : geo.meshes()) {
                    const int orientation = geo.oriented(mesh1,mesh2);
                    if (orientation != 0){
                        // D*_23 or D*_33
                        operatorDstar(mesh2,mesh1,transmat,K*orientation,gauss_order);
                        if (mesh1==mesh2) {
                            // I_33
                            operatorP1P0(mesh1,transmat,-0.5*orientation);
                        } else {
                            // S_23
                            operatorS(mesh2,mesh1,transmat,geo.sigma_inv(mesh1,mesh2)*(-1.0*K*orientation),gauss_order);
                        }
                    }
                }

        for ( size_t ielec = 0; ielec < n_sensors; ++ielec) {
            Triangles tris = electrodes.getInjectionTriangles(ielec);
            for ( Triangles::const_iterator tit = tris.begin(); tit != tris.end(); ++tit) {
                // to ensure exactly no accumulation of currents. w = elec_area/tris_area (~= 1)
                double inv_area = electrodes.getWeights()(ielec);
                // if no radius is given, we assume the user wants to specify an intensity not a density of current
                if ( almost_equal(electrodes.getRadii()(ielec), 0.) ) {
                    inv_area = 1./tit->area();
                }
                for ( size_t i = 0; i < (geo.nb_parameters() - geo.nb_current_barrier_triangles()); ++i) {
                    mat(i, ielec) += transmat(tit->index(), i) * inv_area;
                }
            }
        }
    }

    DipSource2InternalPotMat::DipSource2InternalPotMat(const Geometry& geo,const Matrix& dipoles,const Matrix& points,const std::string& domain_name) {

        Matrix& mat = *this;

        // Points with one more column for the index of the domain they belong

        std::vector<const Domain*> points_domain;
        std::vector<Vect3>   points_;
        for (unsigned i=0; i<points.nlin(); ++i) {
            const Domain& domain = geo.domain(Vect3(points(i,0),points(i,1),points(i,2)));
            if (domain.conductivity()!=0.0) {
                points_domain.push_back(&domain);
                points_.push_back(Vect3(points(i,0),points(i,1),points(i,2)));
            } else {
                std::cerr << " DipSource2InternalPot: Point [ " << points.getlin(i);
                std::cerr << "] is outside the head. Point is dropped." << std::endl;
            }
        }
        const double K = 1.0/(4*Pi);
        mat = Matrix(points_.size(), dipoles.nlin());
        mat.set(0.0);

        for (unsigned iDIP=0; iDIP<dipoles.nlin(); ++iDIP) {
            const Vect3 r0(dipoles(iDIP,0), dipoles(iDIP,1), dipoles(iDIP,2));
            const Vect3  q(dipoles(iDIP,3), dipoles(iDIP,4), dipoles(iDIP,5));

            const Domain& domain = (domain_name=="") ? geo.domain(r0) : geo.domain(domain_name);
            const double  cond   = domain.conductivity();

            static analyticDipPot anaDP;
            anaDP.init(q, r0);
            for (unsigned iPTS=0; iPTS<points_.size(); ++iPTS)
                if (points_domain[iPTS]==&domain)
                    mat(iPTS, iDIP) += K/cond*anaDP.f(points_[iPTS]);
        }
    }
}
