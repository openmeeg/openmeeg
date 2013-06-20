/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
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

#if WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#include "vector.h"
#include "matrix.h"
#include "danielsson.h"
#include "operators.h"
#include "assemble.h"
#include "sensors.h"
#include <fstream>

namespace OpenMEEG {

    using namespace std;

    void assemble_SurfSourceMat(Matrix& mat, const Geometry& geo, const Mesh& sources, const int gauss_order) {
/*
        const int newsize = geo.size() - geo.end()->nb_triangles();
        mat = Matrix(newsize, sources.nb_vertices());
        mat.set(0.0);

        const unsigned nVertexSources    = sources.nb_vertices();
        const unsigned nVertexFirstLayer = geo.begin()->nb_vertices();
        const unsigned nFacesFirstLayer  = geo.begin()->nb_triangles();
        std::cout << std::endl << "assemble SurfSourceMat with " << nVertexSources << " sources" << std::endl << std::endl;

        //  First block is nVertexFistLayer*nVertexSources.

        operatorN(*geo.begin(), sources, mat, 0, 0, gauss_order);

        //  Second block is nFacesFistLayer*nVertexSources.

        operatorD(*geo.begin(), sources, mat, static_cast<int>(nVertexFirstLayer), 0, gauss_order);

        //  First block*=(-1/sigma_inside).

        const double K = 1.0/(4.*M_PI);
        const double s1i = geo.sigma("Brain");

        mult(mat, nVertexFirstLayer, 0, nVertexFirstLayer + nFacesFirstLayer-1, nVertexSources-1, -K/s1i);
        mult(mat, 0, 0, nVertexFirstLayer-1, nVertexSources-1, K); 
        */
    }

    SurfSourceMat::SurfSourceMat(const Geometry& geo, const Mesh& sources, const int gauss_order) {
        assemble_SurfSourceMat(*this, geo, sources, gauss_order);
    }

    void assemble_DipSourceMat(Matrix& rhs, const Geometry& geo, const Matrix& dipoles,
            const int gauss_order, const bool adapt_rhs, const bool dipoles_innermost) {

        const double K         = 1.0/(4.*M_PI);
        const int    newsize   = geo.size()-geo.end()->nb_triangles();
        const size_t n_dipoles = dipoles.nlin();
        rhs = Matrix(newsize, n_dipoles);
        rhs.set(0.);

        //  First block is nVertexFistLayer.

        Vector rhs_col(rhs.nlin());
        for (size_t s = 0; s < n_dipoles; s++) {
            PROGRESSBAR(s, n_dipoles);
            const Vect3 r(dipoles(s, 0), dipoles(s, 1), dipoles(s, 2));
            const Vect3 q(dipoles(s, 3), dipoles(s, 4), dipoles(s, 5));

            const Domain domain = (dipoles_innermost && geo.domains()[0].innermost()) ? geo.domains()[0] : geo.get_domain(r); // TODO check
            // TODO add conditions on Air domain : dipoles inside scalp ?
            const double sigma  = domain.sigma();

            rhs_col.set(0.);
            // iterate over the domain's interfaces (half-spaces)
            for (Domain::const_iterator hit = domain.begin(); hit != domain.end(); hit++ ) {
                // iterate over the meshes of the interface
                for (Interface::const_iterator mit = hit->interface().begin(); mit != hit->interface().end(); mit++ ) {
                    //  Treat the mesh.
                    operatorDipolePotDer(r, q, **mit, rhs_col, gauss_order, adapt_rhs);

                    for (Mesh::const_vertex_iterator vit = (*mit)->vertex_begin(); vit!= (*mit)->vertex_end(); mit++) {
                        rhs_col((*vit)->index()) *= (hit->inside())?K:-K; // TODO check signs
                    }

                    if (!(*mit)->outermost()) {
                        operatorDipolePot(r, q, **mit, rhs_col, gauss_order, adapt_rhs);

                        for (Mesh::const_iterator tit = (*mit)->begin(); tit!= (*mit)->end(); tit++) {
                            rhs_col(tit->index()) *= (hit->inside())?-K/sigma:(K/sigma);
                        }
                    }
                }
            }
            rhs.setcol(s, rhs_col);
        }
    }

    DipSourceMat::DipSourceMat(const Geometry& geo, const Matrix& dipoles, const int gauss_order,
                                const bool adapt_rhs, const bool dipoles_innermost) {

        assemble_DipSourceMat(*this, geo, dipoles, gauss_order, adapt_rhs, dipoles_innermost);
    }

    void assemble_EITSourceMat(Matrix& mat, const Geometry& geo, const Matrix& positions, const int gauss_order) {

        /*
        //  A Matrix to be applied to the scalp-injected current to obtain the Source Term of the EIT foward problem.

        const double K = 1.0/(4.*M_PI);
        const int newsize = geo.size()-geo.end()->nb_triangles();
        mat = Matrix(newsize, positions.nlin());

        //  transmat = a big SymMatrix of which mat = part of its transpose.
        SymMatrix transmat(geo.size());
        transmat.set(0.0);

        int offset0, offset1, offset2, offset3, offset4;

        // We iterate over the meshes (or pair of domains)
        for (Geometry::const_iterator mit = geo.begin(); mit != geo.end(); mit++) {
            offset0 = offset;
            offset1 = offset0 + mit->nb_vertices();
            offset2 = offset1 + mit->nb_triangles();
            offset3 = offset2 + geo.getM(c+1).nb_vertices();  // TODO
            offset4 = offset3 + geo.getM(c+1).nb_triangles(); // TODO
            offset  = offset2;
        }

        //  Compute S.

        operatorS(geo.getM(c+1), geo.getM(c), transmat, offset3, offset1, gauss_order);
        mult(transmat, offset3, offset1, offset4, offset2, K/geo.sigma_in(c+1));

        //  First compute D, then it will be transposed.

        operatorD(geo.getM(c+1), geo.getM(c), transmat, offset3, offset0, gauss_order);
        mult(transmat, offset3, offset0, offset4, offset1, -K);
        operatorD(geo.getM(c+1), geo.getM(c+1), transmat, offset3, offset2, gauss_order);
        mult(transmat, offset3, offset2, offset4, offset3, -2.0*K);
        operatorP1P0(geo.getM(c+1), transmat, offset3, offset2);
        mult(transmat, offset3, offset2, offset4, offset3, -0.5);

        //  Extracting the transpose of the last block of lines of transmat.
        //  Transposing the Matrix.

        for (unsigned ielec = 0; ielec < positions.nlin(); ++ielec) {
            const Vect3 current_position(positions(ielec, 0), positions(ielec, 1), positions(ielec, 2));
            Vect3 current_alphas; //not used here
            int current_nearest_triangle; //    To hold the index of the closest triangle to electrode.
            dist_point_mesh(current_position, geo.getM(geo.nb()-1), current_alphas, current_nearest_triangle);
            const double inv_area = 1.0/geo.getM(geo.nb()-1).triangle(current_nearest_triangle).area();
            for (int i = 0; i < newsize; ++i) {
                mat(i, ielec) = transmat(newsize+current_nearest_triangle, i)*inv_area;
            }
        }
        */
    }

    EITSourceMat::EITSourceMat(const Geometry& geo, Sensors& electrodes, const int gauss_order) {
        assemble_EITSourceMat(*this, geo, electrodes.getPositions(), gauss_order);
    }

    void assemble_DipSource2InternalPotMat(Matrix& mat, const Geometry& geo, const Matrix& dipoles,
                                           const Matrix& points, const bool dipoles_innermost)     {
        // Points with one more column for the index of the domain they belong
        std::vector<Domain> points_domain;
        for (unsigned i = 0; i < points.nlin(); ++i) {
            points_domain.push_back(geo.get_domain(Vect3(points(i, 0), points(i, 1), points(i, 2))));
        }
        const double K = 1.0/(4.*M_PI);
        mat = Matrix(points.nlin(), dipoles.nlin());
        mat.set(0.0);

        for (unsigned iDIP = 0; iDIP < dipoles.nlin(); ++iDIP) {
            const Vect3 r0(dipoles(iDIP, 0), dipoles(iDIP, 1), dipoles(iDIP, 2));
            const Vect3  q(dipoles(iDIP, 3), dipoles(iDIP, 4), dipoles(iDIP, 5));

            const Domain domain = (dipoles_innermost && geo.domains()[0].innermost()) ? geo.domains()[0] : geo.get_domain(r0); // TODO check
            const double sigma  = domain.sigma();

            static analyticDipPot anaDP;
            anaDP.init(q, r0);
            for (size_t iPTS = 0; iPTS < points.nlin(); ++iPTS) {
                if (points_domain[iPTS] == domain) {
                    const Vect3 r(points(iPTS, 0), points(iPTS, 1), points(iPTS, 2));
                    mat(iPTS, iDIP) = K/sigma*anaDP.f(r);
                }
            }
        }
    }

    DipSource2InternalPotMat::DipSource2InternalPotMat(const Geometry& geo, const Matrix& dipoles,
                                                       const Matrix& points, const bool dipoles_innermost)
    {
        assemble_DipSource2InternalPotMat(*this, geo, dipoles, points, dipoles_innermost);
    }
}
