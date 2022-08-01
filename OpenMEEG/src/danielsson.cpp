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

#include <danielsson.h>
/* implements an algorithm  proposed in Danielsson, P.-E. Euclidean Distance Mapping. Computer Graphics and Image
Processing 14, 3 (Nov. 1980), 227-248.
*/
namespace OpenMEEG
{

    // Distance from p to a triangle (3 pointers to vertex)
    // alpha-> barycentric coordinates of closest point
    // sum(alpha_i)=1
    // inside: closest point is inside (alpha_i!=0 for all i)

    // Auxilary Fn : nb vertices left (for the others alpha=0)

    using namespace std;

    static double dpc(const Vect3& p,const Triangle& triangle,Vect3& alphas,unsigned nb,int* idx,bool& inside) {
        if (nb==1) {
            alphas(idx[0]) = 1.0;
            return (p-triangle.vertex(idx[0])).norm();
        }
        // Solves H=sum(alpha_i A_i), sum(alpha_i)=1, et HM.(A_i-A_0)=0
        Vect3 A0Ai[3]; // A_i-A_0
        for (unsigned i=1; i<nb; ++i)
            A0Ai[i] = triangle.vertex(idx[i])-triangle.vertex(idx[0]);

        Vect3 A0M = p-triangle.vertex(idx[0]); // M-A_0
        if (nb==2) {
            alphas(idx[1]) = dotprod(A0M,A0Ai[1])/dotprod(A0Ai[1],A0Ai[1]);
            alphas(idx[0]) = 1.0-alphas(idx[1]);
        } else if (nb==3) {
            // direct inversion (2x2 linear system)
            const double a00 = dotprod(A0Ai[1],A0Ai[1]);
            const double a10 = dotprod(A0Ai[1],A0Ai[2]);
            const double a11 = dotprod(A0Ai[2],A0Ai[2]);
            const double b0 = dotprod(A0M,A0Ai[1]);
            const double b1 = dotprod(A0M,A0Ai[2]);
            const double d = a00*a11-a10*a10;
            om_error(d!=0);
            alphas(idx[1]) = (b0*a11-b1*a10)/d;
            alphas(idx[2]) = (a00*b1-a10*b0)/d;
            alphas(idx[0]) = 1.0-alphas(idx[1])-alphas(idx[2]);
        } else {
            // 3 unknowns or more -> solve system
            //  Ax=b with: A(i, j)=A0Ai.AjA0, x=(alpha_1, alpha_2, ...), b=A0M.A0Ai
            cerr << "Error : dim>=4 in danielsson !" << endl;
            exit(0);
        }
        // If alpha_i<0 -> brought to 0 and recursion
        // NB: also takes care of alpha > 1 because if alpha_i>1 then alpha_j<0 for at least one j
        for (unsigned i=0; i<nb; ++i) {
            if (alphas(idx[i])<0) {
                inside = false;
                alphas(idx[i]) = 0;
                swap(idx[i], idx[nb-1]);
                return dpc(p, triangle, alphas, nb-1, idx, inside);
            }
        }
        // Sinon: distance HM
        Vect3 MH = -A0M;
        for (unsigned i=1; i<nb; ++i)
            MH += alphas(idx[i]) * A0Ai[i];
        return MH.norm();
    }

    // Main Function
    double dist_point_triangle(const Vect3& p,const Triangle& triangle,Vect3& alphas,bool& inside) {
        int idx[3] = { 0, 1, 2 };
        inside = true;
        return dpc(p,triangle,alphas,3,idx,inside);
    }

    static inline int sgn(double s) {
        return (s>0) ? 1 : (s<0) ? -1 : 0;
    }

    std::tuple<double,const Triangle&,const Mesh&>
    dist_point_interface(const Vect3& p,const Interface& interface,Vect3& alphas) {
        double distmin = std::numeric_limits<double>::max();
        const Mesh*     nearest_mesh     = nullptr;
        const Triangle* nearest_triangle = nullptr;
        for (const auto& omesh : interface.oriented_meshes())
            for (const auto& triangle : omesh.mesh().triangles()) {
                bool  inside;
                Vect3 alphasLoop;
                const double distance = dist_point_triangle(p,triangle,alphasLoop,inside);
                if (distance<distmin) {
                    distmin = distance;
                    alphas  = alphasLoop;
                    nearest_triangle = &triangle;
                    nearest_mesh = &omesh.mesh();
                }
            }
        return { distmin, *nearest_triangle, *nearest_mesh };;
    }

    // Find the closest triangle on the interfaces that touches 0 conductivity

    std::tuple<double,const Triangle&,const Mesh&,const Interface&>
    dist_point_geom(const Vect3& p,const Geometry& g,Vect3& alphas) {

        double distmin = std::numeric_limits<double>::max();
        const Mesh*      nearest_mesh      = nullptr;
        const Interface* nearest_interface = nullptr;
        const Triangle*  nearest_triangle  = nullptr;

        for(const auto& domain : g.domains())
            if (domain.conductivity()==0.0)
                for (const auto& boundary : domain.boundaries()) {
                    const auto& res = dist_point_interface(p,boundary.interface(),alphas);
                    const double distance = std::get<0>(res);
                    if (distance<distmin) {
                        nearest_interface = &boundary.interface();
                        nearest_mesh = &std::get<2>(res);
                        distmin = distance;
                        nearest_triangle = &std::get<1>(res);
                    }
                }

        return { distmin, *nearest_triangle, *nearest_mesh, *nearest_interface };
    }
}
