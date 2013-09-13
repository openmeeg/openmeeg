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
    void deflat(T& M, const Interface& i, double coef) 
    {
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

    void assemble_HM(const Geometry& geo, SymMatrix& mat, const unsigned gauss_order) 
    {
        mat = SymMatrix((geo.size()-geo.outermost_interface().nb_triangles()));
        mat.set(0.0);
        double K = 1.0 / (4.0 * M_PI);

        // We iterate over the meshes (or pair of domains) to fill the lower half of the HeadMat (since its symmetry)
        for ( Geometry::const_iterator mit1 = geo.begin(); mit1 != geo.end(); ++mit1) {

            for ( Geometry::const_iterator mit2 = geo.begin(); (mit2 != (mit1+1)); ++mit2) {

                // if mit1 and mit2 communicate, i.e they are used for the definition of a common domain
                const double orientation = geo.oriented(*mit1, *mit2); // equals  0, if they don't have any domains in common
                                                                       // equals  1, if they are both oriented toward the same domain
                                                                       // equals -1, if they are not

                if ( std::abs(orientation) > 10.*std::numeric_limits<double>::epsilon() ) {

                    double Scoeff =   orientation * geo.sigma_inv(*mit1, *mit2) * K;
                    double Dcoeff = - orientation * geo.indicator(*mit1, *mit2) * K;
                    double Ncoeff;

                    if ( !(mit1->outermost() || mit2->outermost()) ) {
                        // Computing S block first because it's needed for the corresponding N block
                        operatorS(*mit1, *mit2, mat, Scoeff, gauss_order);
                        Ncoeff = geo.sigma(*mit1, *mit2)/geo.sigma_inv(*mit1, *mit2);
                    } else {
                        Ncoeff = orientation * geo.sigma(*mit1, *mit2) * K;
                    }

                    if ( !mit1->outermost() ) {
                        // Computing D block
                        operatorD(*mit1, *mit2, mat, Dcoeff, gauss_order);
                    }
                    if ( ( *mit1 != *mit2 ) && ( !mit2->outermost() ) ) {
                        // Computing D* block
                        operatorD(*mit1, *mit2, mat, Dcoeff, gauss_order, true);
                    }

                    // Computing N block
                    operatorN(*mit1, *mit2, mat, Ncoeff, gauss_order);
                }
            }
        }

        // Deflate the last diagonal block of 'mat' : (in order to have a zero-mean potential for the outermost interface)
        const Interface i = geo.outermost_interface();
        unsigned i_first = (*i.begin()->mesh().vertex_begin())->index();
        deflat(mat, i, mat(i_first, i_first) / (geo.outermost_interface().nb_vertices()));
    }

    void assemble_Surf2Vol(const Geometry& geo, Matrix& mat, const std::map<const Domain, Vertices> m_points) 
    {
        const double K = 1.0/(4.0*M_PI);

        unsigned size = 0; // total number of inside points
        for ( std::map<const Domain, Vertices>::const_iterator dvit = m_points.begin(); dvit != m_points.end(); ++dvit) {
            size += dvit->second.size();
        }

        mat = Matrix(size, (geo.size() - geo.outermost_interface().nb_triangles()));
        mat.set(0.0);

        for ( std::map<const Domain, Vertices>::const_iterator dvit = m_points.begin(); dvit != m_points.end(); ++dvit) {
            for ( Geometry::const_iterator mit = geo.begin(); mit != geo.end(); ++mit) {
                int orientation = dvit->first.mesh_orientation(*mit);
                if ( orientation != 0 ) {
                    operatorDinternal(*mit, mat, dvit->second, orientation * -1. * K);
                    if ( !mit->outermost() ) {
                        operatorSinternal(*mit, mat, dvit->second, orientation * K / geo.sigma(dvit->first));
                    }
                }
            }
        }
    }

    HeadMat::HeadMat(const Geometry& geo, const unsigned gauss_order)
    {
        assemble_HM(geo, *this, gauss_order);
    }

    Surf2VolMat::Surf2VolMat(const Geometry& geo, const Matrix& points) 
    {
        std::map<const Domain, Vertices> m_points;

        unsigned index = 0;
        // Find the points per domain and generate the indices for the m_points
        for ( unsigned i = 0; i < points.nlin(); ++i) {
            const Domain domain = geo.domain(Vect3(points(i, 0), points(i, 1), points(i, 2)));
            if ( domain.name() == "Air" ) {
                std::cerr << " Surf2Vol: Point [ " << points.getlin(i);
                std::cerr << "] is outside the head. Point is dropped." << std::endl;
            } else {
                m_points[domain].push_back(Vertex(points(i, 0), points(i, 1), points(i, 2), index++));
            }
        }

        assemble_Surf2Vol(geo, *this, m_points);
    }

} // namespace OpenMEEG
