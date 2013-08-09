/* OpenMEEG

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

#include <matrix.h>
#include <symmatrix.h>
#include <geometry.h>
#include <operators.h>
#include <assemble.h>

namespace OpenMEEG {

    template<class T>
    void deflat(T& M, const Interface& i, double coef) {
        // deflate the Matrix
        for ( Interface::const_iterator omit = i.begin(); omit != i.end(); ++omit) {
            for ( Mesh::const_vertex_iterator vit1 = omit->mesh().vertex_begin(); vit1 != omit->mesh().vertex_end(); ++vit1) {
                #pragma omp parallel for
                for ( Mesh::const_vertex_iterator vit2 = vit1; vit2 < omit->mesh().vertex_end(); ++vit2) {
                    M((*vit1)->index(), (*vit2)->index()) += coef;
                }
            }
        }
    }

    void assemble_HM(const Geometry& geo, SymMatrix& mat, const unsigned gauss_order) {

        mat = SymMatrix((geo.size()-geo.outermost_interface().nb_triangles()));
        mat.set(0.0);
        double K = 1.0 / (4.0 * M_PI);

        // We iterate over the meshes (or pair of domains) to fill the lower half of the headmat (since its symmetry)
        for ( Geometry::const_iterator mit1 = geo.begin(); mit1 != geo.end(); ++mit1) {

            for ( Geometry::const_iterator mit2 = geo.begin(); (mit2 != (mit1+1)); ++mit2) {

                // if mit1 and mit2 communicate, i.e they are used for the definition of a domain
                const double orientation = geo.oriented(*mit1, *mit2); // equals  0, if they don't have any domains in common
                                                                       // equals  1, if they are both oriented toward the same domain
                                                                       // equals -1, if they are not

                if ( std::abs(orientation) > 10.*std::numeric_limits<double>::epsilon() ) {

                    if ( !(mit1->outermost() || mit2->outermost()) ) {
                        // Computing S block first because it's needed for the corresponding N block
                        double Scoeff =   orientation * geo.sigma_inv(*mit1, *mit2) * K;
                        operatorS(*mit1, *mit2, mat, Scoeff, gauss_order);
                    }
                    if ( !mit1->outermost() ) {
                        // Computing D block
                        double Dcoeff = - orientation * geo.indicatrice(*mit1, *mit2) * K;
                        operatorD(*mit1, *mit2, mat, Dcoeff, gauss_order);
                    }
                    if ( ( *mit1 != *mit2 ) && ( !mit2->outermost() ) ) {
                        // Computing D* block
                        double Dcoeff = - orientation * K;
                        operatorD(*mit1, *mit2, mat, Dcoeff, gauss_order, true);
                    }

                    // Computing N block
                    double Ncoeff = orientation * geo.sigma(*mit1, *mit2) * K;
                    operatorN(*mit1, *mit2, mat, Ncoeff, gauss_order);
                }
            }
        }

        // Deflate the last diagonal block of 'mat' : (in order to have a zero-mean potential for the outermost interface)
        const Interface i = geo.outermost_interface();
        unsigned i_first = (*i.begin()->mesh().vertex_begin())->index();
        deflat(mat, i, mat(i_first, i_first) / (geo.outermost_interface().nb_vertices()));
    }

    void assemble_Surf2Vol(const unsigned N, const Geometry& geo, Matrix& mat, const std::vector<Matrix>& points) 
    {
        const double K = 1.0/(4.0*M_PI);

        mat = Matrix(N, (geo.size()-geo.outermost_interface().nb_triangles()));
        mat.set(0.0);

// TODO
#if 0
        // We iterate over the meshes (or pair of domains) to fill the lower half of the headmat (since its symmetry)
        for (Geometry::const_iterator mit1 = geo.begin(); mit1 != geo.end(); ++mit1) {

            for (Geometry::const_iterator mit2 = geo.begin(); (mit2 != (mit1+1)); ++mit2) {

                // if mit1 and mit2 communicate, i.e they are used for the definition of a domain
                const double orientation = geo.oriented(*mit1, *mit2); // equals  0, if they don't have any domains in common

        unsigned offset = 0;
        unsigned offsetA0 = 0;
        unsigned c = 0;
        for (Geometry::const_iterator mit = geo.begin(); mit != geo.end(); ++mit, ++c) {
            const unsigned offset0 = offset;
            const unsigned offsetA = offsetA0;
            const unsigned offsetB = offsetA + points[c].nlin();
            const unsigned offset1 = offset  + mit->nb_vertices();
            const unsigned offset2 = offset1 + mit->size();
            const unsigned offset3 = offset2 + (mit+1)->nb_vertices();
            const unsigned offset4 = offset3 + (mit+1)->size();

            // compute DI, i block if necessary. // TODO check les orientations
            if ( c==0) {
                operatorDinternal(*mit, mat, offsetA, offset0, points[c]);
                mult(mat, offsetA, offset0, offsetB, offset1, -K);
            }

            // compute DI+1, i block.
            operatorDinternal(*mit, mat, offsetB, offset0, points[c+1]);
            mult(mat, offsetB, offset0, offsetB+points[c+1].nlin(), offset1, K);

            // compute DI+1, i+1 block
            operatorDinternal(*mit, mat, offsetB, offset2, points[c+1]);
            mult(mat, offsetB, offset2, offsetB+points[c+1].nlin(), offset3,-K);

            // compute SI, i block if necessary.
            if ( c == 0 ) {
                operatorSinternal(*mit, mat, offsetA, offset1, points[c]);
                mult(mat, offsetA, offset1, offsetA+points[c].nlin(), offset2, K/geo.sigma_in(c));
            }

            // compute SI+1, i block
            const double inv_sig = K/geo.sigma_in(c+1);
            operatorSinternal(*mit, mat, offsetB, offset1, points[c+1]);
            mult(mat, offsetA+points[c].nlin(), offset1, offsetB+points[c+1].nlin(), offset2, -inv_sig);

            if ( c < geo.nb_domains()-2 ) {
                // compute SI+1, i block
                operatorSinternal(*(mit+1), mat, offsetB, offset3, points[c+1]);
                mult(mat, offsetB, offset3, offsetB+points[c+1].nlin(), offset4, inv_sig);
            }
            offset   = offset2;
            offsetA0 = offsetA + points[c].nlin();
        }
#endif
    }

    HeadMat::HeadMat (const Geometry& geo, const unsigned gauss_order)
    {
        assemble_HM(geo, *this, gauss_order);
    }

    Surf2VolMat::Surf2VolMat(const Geometry& geo, const Matrix& points) 
    {
#if 0
        //  Find and count the points per domain.
        std::vector<unsigned> labels(points.nlin());
        std::vector<unsigned> nb_pts_per_dom(geo.nb_domains(), 0);
        for ( unsigned i = 0; i < points.nlin(); ++i) {
            unsigned domain = geo.domain(Vect3(points(i, 0), points(i, 1), points(i, 2)));
            if ( domain >= geo.nb_domains() ) {
                std::cerr << " Surf2Vol: Point " << points.getlin(i);
                std::cerr << " is outside the head. Point is considered to be in the scalp." << std::endl;
                domain = geo.nb_domains()-1;
            }
            labels[i] = domain;
            ++nb_pts_per_dom[domain];
        }

        //  Split the point array into multiple arrays (one array per domain).

        std::vector<Matrix> vect_PtsInDom(geo.nb_domains());
        for ( unsigned c = 0; c < static_cast<unsigned>(geo.nb_domains()); ++c) {
            vect_PtsInDom[c] = Matrix(nb_pts_per_dom[c], 3);
            for ( unsigned ipt = 0, iptd = 0; ipt < points.nlin(); ++ipt) { // get the points in the domain c
                if ( labels[ipt] == c ) {
                    for ( unsigned ic = 0; ic < 3; ++ic) {
                        vect_PtsInDom[c](iptd, ic) = points(ipt, ic);
                    }
                    ++iptd;
                }
            }
        }

        assemble_Surf2Vol(points.nlin(), geo, *this, vect_PtsInDom);
#endif
    }

} // namespace OpenMEEG
