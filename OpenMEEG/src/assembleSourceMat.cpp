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

#if WIN32
#define _USE_MATH_DEFINES
#endif
#include <cmath>

#include <vector.h>
#include <matrix.h>
#include <danielsson.h>
#include <operators.h>
#include <assemble.h>
#include <sensors.h>
#include <fstream>

namespace OpenMEEG {

    void assemble_SurfSourceMat(Matrix& mat, const Geometry& geo, Mesh& mesh_source, const unsigned gauss_order) 
    {
        mat = Matrix((geo.size()-geo.nb_current_barrier_triangles()), mesh_source.nb_vertices());
        mat.set(0.0);

        // check if no overlapping between the geometry and the source mesh
        bool OK = geo.check(mesh_source);
        if ( !OK ) {
            std::cerr << "Error: source mesh overlapps the geometry" << std::endl;
            return;
        } // then the mesh is included in a domain of the geometry

        const Domain d     = geo.domain(**mesh_source.vertex_begin()); 
        const double sigma = d.sigma();
        const double K     = 1.0/(4.*M_PI);
        
        // We here set it as an outermost (to tell _operarorN it doesn't belong to the geometry)
        mesh_source.outermost() = true;
        mesh_source.current_barrier() = true;

        std::cout << std::endl << "assemble SurfSourceMat with " << mesh_source.nb_vertices() << " mesh_source located in domain \"" << d.name() << "\"." << std::endl << std::endl;

        for ( Domain::const_iterator hit = d.begin(); hit != d.end(); ++hit) {
            for ( Interface::const_iterator omit = hit->interface().begin(); omit != hit->interface().end(); ++omit) {
                // First block is nVertexFistLayer*mesh_source.nb_vertices()
                double coeffN = (hit->inside())?K * omit->orientation() : omit->orientation() * -K;
                operatorN( omit->mesh(), mesh_source, mat, coeffN, gauss_order);
                // Second block is nFacesFistLayer*mesh_source.nb_vertices()
                double coeffD = (hit->inside())?-omit->orientation() * K / sigma : omit->orientation() * K / sigma;
                operatorD(omit->mesh(), mesh_source, mat, coeffD, gauss_order,false);
            }
        }
    }

    SurfSourceMat::SurfSourceMat(const Geometry& geo, Mesh& mesh_source, const unsigned gauss_order) 
    {
        assemble_SurfSourceMat(*this, geo, mesh_source, gauss_order);
    }

    void assemble_DipSourceMat(Matrix& rhs, const Geometry& geo, const Matrix& dipoles,
            const unsigned gauss_order, const bool adapt_rhs, const std::string& domain_name = "") 
    {
        const size_t size      = geo.size()-geo.nb_current_barrier_triangles();
        const size_t n_dipoles = dipoles.nlin();

        rhs = Matrix(size,n_dipoles);
        rhs.set(0.0);

        Vector rhs_col(rhs.nlin());
        for (unsigned s=0; s<n_dipoles; ++s) {
            PROGRESSBAR(s,n_dipoles);
            const Vect3 r(dipoles(s,0),dipoles(s,1),dipoles(s,2));
            const Vect3 q(dipoles(s,3),dipoles(s,4),dipoles(s,5));

            const Domain domain = (domain_name=="") ? geo.domain(r) : geo.domain(domain_name);

            //  Only consider dipoles in non-zero conductivity domain.

            const double sigma = domain.sigma();
            if (sigma!=0.0) {
                rhs_col.set(0.0);
                const double K = 1.0/(4.*M_PI);
                //  Iterate over the domain's interfaces (half-spaces)
                for (Domain::const_iterator hit=domain.begin(); hit!=domain.end(); ++hit) {
                    //  Iterate over the meshes of the interface
                    for (Interface::const_iterator omit=hit->interface().begin(); omit!=hit->interface().end(); ++omit) {
                        //  Treat the mesh.
                        const double coeffD = ((hit->inside()) ? K : -K)*omit->orientation();
                        operatorDipolePotDer(r,q,omit->mesh(),rhs_col,coeffD,gauss_order,adapt_rhs);

                        if (!omit->mesh().current_barrier()) {
                            const double coeff = -coeffD/sigma;;
                            operatorDipolePot(r,q,omit->mesh(),rhs_col,coeff,gauss_order,adapt_rhs);
                        }
                    }
                }
                rhs.setcol(s,rhs_col);
            }
        }
    }

    DipSourceMat::DipSourceMat(const Geometry& geo, const Matrix& dipoles, const unsigned gauss_order,
                               const bool adapt_rhs, const std::string& domain_name)
    {
        assemble_DipSourceMat(*this, geo, dipoles, gauss_order, adapt_rhs, domain_name);
    }

    void assemble_EITSourceMat(Matrix& mat, const Geometry& geo, const Sensors& electrodes, const unsigned gauss_order)
    {
        //  A Matrix to be applied to the scalp-injected current to obtain the Source Term of the EIT foward problem.
        // following article BOUNDARY ELEMENT FORMULATION FOR ELECTRICAL IMPEDANCE TOMOGRAPHY
        // (eq.14 (do not look at eq.16 since there is a mistake: D_23 -> S_23))
        // rhs = [0 ... 0  -D*_23  sigma_3^(-1)S_23  -I_33/2.+D*_33]

        size_t n_sensors = electrodes.getNumberOfSensors();

        const double K = 1.0/(4.*M_PI);

        //  transmat = a big SymMatrix of which mat = part of its transpose.
        SymMatrix transmat(geo.size());
        transmat.set(0.0);
        mat = Matrix((geo.size()-geo.nb_current_barrier_triangles()), n_sensors);
        mat.set(0.0);

        for (Geometry::const_iterator mit0 = geo.begin(); mit0 != geo.end(); ++mit0) {
            if (mit0->current_barrier()) {
                for (Geometry::const_iterator mit1 = geo.begin(); mit1 != geo.end(); ++mit1) {
                    const int orientation = geo.oriented(*mit0,*mit1);
                    if (orientation != 0){
                        // D*_23 or D*_33
                        operatorD(*mit1, *mit0, transmat, K*orientation, gauss_order, true);
                        if (*mit0==*mit1) {
                            // I_33
                            operatorP1P0(*mit0, transmat, -0.5*orientation);
                        } else {
                            // S_23
                            operatorS(*mit1, *mit0, transmat, geo.sigma_inv(*mit0,*mit1)*(-1.0*K*orientation), gauss_order);
                        }
                    }
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
                for ( size_t i = 0; i < (geo.size() - geo.nb_current_barrier_triangles()); ++i) {
                    mat(i, ielec) += transmat(tit->index(), i) * inv_area;
                }
            }
        }
    }

    EITSourceMat::EITSourceMat(const Geometry& geo, const Sensors& electrodes, const unsigned gauss_order) 
    {
        assemble_EITSourceMat(*this, geo, electrodes, gauss_order);
    }

    void assemble_DipSource2InternalPotMat(Matrix& mat, const Geometry& geo, const Matrix& dipoles,
                                           const Matrix& points, const std::string& domain_name)     
    {
        // Points with one more column for the index of the domain they belong
        std::vector<Domain> points_domain;
        std::vector<Vect3>  points_;
        for ( unsigned i = 0; i < points.nlin(); ++i) {
            const Domain& d = geo.domain(Vect3(points(i, 0), points(i, 1), points(i, 2)));
            if ( d.sigma() != 0.0 ) {
                points_domain.push_back(d);
                points_.push_back(Vect3(points(i, 0), points(i, 1), points(i, 2)));
            }
            else {
                std::cerr << " DipSource2InternalPot: Point [ " << points.getlin(i);
                std::cerr << "] is outside the head. Point is dropped." << std::endl;
            }
        }
        const double K = 1.0/(4.*M_PI);
        mat = Matrix(points_.size(), dipoles.nlin());
        mat.set(0.0);

        for ( unsigned iDIP = 0; iDIP < dipoles.nlin(); ++iDIP) {
            const Vect3 r0(dipoles(iDIP, 0), dipoles(iDIP, 1), dipoles(iDIP, 2));
            const Vect3  q(dipoles(iDIP, 3), dipoles(iDIP, 4), dipoles(iDIP, 5));

            Domain domain;
            if ( domain_name == "" ) {
                domain = geo.domain(r0);
            } else {
                domain = geo.domain(domain_name);
            }
            const double sigma  = domain.sigma();

            static analyticDipPot anaDP;
            anaDP.init(q, r0);
            for ( unsigned iPTS = 0; iPTS < points_.size(); ++iPTS) {
                if ( points_domain[iPTS] == domain ) {
                    mat(iPTS, iDIP) += K/sigma*anaDP.f(points_[iPTS]);
                }
            }
        }
    }

    DipSource2InternalPotMat::DipSource2InternalPotMat(const Geometry& geo, const Matrix& dipoles,
                                                       const Matrix& points, const std::string& domain_name)
    {
        assemble_DipSource2InternalPotMat(*this, geo, dipoles, points, domain_name);
    }
}
