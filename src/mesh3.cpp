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

#include <OpenMEEGConfigure.h>

#if WIN32
#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_DEPRECATE 1
#endif
#include <math.h>
#include <assert.h>
#include "vect3.h"
#include "triangle.h"
#include "mesh3.h"
#include "om_utils.h"
#include "Triangle_triangle_intersection.h"
#include "IOUtils.H"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>


#ifdef WIN32
inline double log2( double n )
{
    // log(n)/log(2) is log2.
    return log( n ) / log( 2.0 );
}
#endif

namespace OpenMEEG {

    Mesh::Mesh(): npts(0) {};

    Mesh::Mesh(const int a, const int b): npts(a), ntrgs(b), pts(new Vect3[npts]),
                            trgs(new Triangle[ntrgs]), links(new intSet[npts]),
                            normals(new Vect3[npts]) { }

    Mesh::Mesh(const Mesh& M): npts(0) { *this = M; }

    Mesh& Mesh::operator= (const Mesh& M) {
        if (this != &M)
            copy(M);
        return *this;
    }

    void Mesh::copy(const Mesh& M) {
        destroy();
        if (M.npts != 0) {
            npts    = M.npts;
            ntrgs   = M.ntrgs;
            pts     = new Vect3[npts];
            trgs    = new Triangle[ntrgs];
            links   = new intSet[npts];
            normals = new Vect3[npts];
            for (int i=0; i < npts; i++) {
                pts[i]     = M.pts[i];
                links[i]   = M.links[i];
                normals[i] = M.normals[i];
            }
            for (int i=0; i < ntrgs; i++)
                trgs[i] = M.trgs[i];
        }
    }

    void Mesh::destroy() {
        // if (npts != 0) {
            // delete [] pts;
            // delete [] trgs;
            // delete [] links;
            // delete [] normals;
// 
            // npts = 0;
            // ntrgs = 0;
        // }
    }

    void Mesh::make_links() {
        for (int i=0; i<npts; i++)
            links[i].clear();

        for (int i=0; i<ntrgs; i++) {
            links[trgs[i].s1()].insert(i);
            links[trgs[i].s2()].insert(i);
            links[trgs[i].s3()].insert(i);
        }
    }

    void Mesh::flip_faces() {
        for( int i = 0; i < ntrgs; ++i )
        {
            Triangle& t = trgs[i];
            int tmp = t[1];
            t[1] = t[0];
            t[0] = tmp;
        }
    }



    /**
     * Append another Mesh to on instance of Mesh
    **/
    void Mesh::append(const Mesh* m) {

        int old_npts = npts;
        int old_ntrgs = ntrgs;

        int new_npts = old_npts + m->nbPts();
        int new_ntrgs = old_ntrgs + m->nbTrgs();

        Vect3 *newPts = new Vect3[new_npts];
        Vect3 *newNormals = new Vect3[new_npts];
        Triangle *newTrgs = new Triangle[new_ntrgs];

        for(int i = 0; i < old_npts; ++i) {
            newPts[i] = point(i);
            newNormals[i] = normal(i);
        }

        for(int i = 0; i < m->nbPts(); ++i) {
            newPts[i + old_npts] = m->point(i);
            newNormals[i + old_npts] = m->normal(i);
        }

        for(int i = 0; i < old_ntrgs; ++i)
            newTrgs[i] = trgs[i];

        for(int i = 0; i < m->nbTrgs(); ++i) {
            const Triangle t = m->triangle(i);
            Triangle* new_t = new Triangle(t[0] + old_npts, t[1] + old_npts, t[2] + old_npts, t.normal());
            newTrgs[i + old_ntrgs] = *new_t;
        }

        destroy(); // Clean old Mesh
        npts = new_npts;
        ntrgs = new_ntrgs;
        pts = newPts;
        trgs = newTrgs;
        normals = newNormals;
        links = new intSet[npts];
        make_links(); // To keep a valid Mesh
        return;
    }

    /**
     * Get the neighboring triangle to triangle (a, b, c) containing edge (a, b)
    **/
    int Mesh::neighbor_triangle(int a, int b, int c) const {
        intSet possible_triangles = links[a];
        intSet::iterator it;
        for(it = possible_triangles.begin(); it != possible_triangles.end(); ++it) {
            Triangle t = triangle(*it);
            if (t.contains(b) && !t.contains(c)) {
                return *it;
            }
        }
        return -1; // Impossible to find neighboring triangle
    }

    /**
     * Print informations about the mesh
    **/
    void Mesh::info() const {

        std::cout << "Mesh Info : " << std::endl;
        std::cout << "\t# points : " << npts << std::endl;
        std::cout << "\t# triangles : " << ntrgs << std::endl;
        std::cout << "\tEuler characteristic : " << npts - 3*ntrgs/2 + ntrgs << std::endl;

        double min_area = trgs[0].area();
        double max_area = trgs[0].area();
        for(int i = 0; i < ntrgs; ++i) {
            min_area = std::min(trgs[i].area(), min_area);
            max_area = std::max(trgs[i].area(), max_area);
        }

        std::cout << "\tMin Area : " << min_area << std::endl;
        std::cout << "\tMax Area : " << max_area << std::endl;
    }

    /**
     * Surface Gradient
    **/
    SparseMatrix Mesh::gradient() const { // TODO a quoi ça sert ?

        SparseMatrix A(3*ntrgs, npts); // nb edges x points
        // loop on triangles
        for (int t=0; t<ntrgs; t++) {
            const Triangle& trg = triangle(t);
            Vect3 points[3] = {point(trg[0]), point(trg[1]), point(trg[2])};
            for(int j=0; j<3; j++) {
                Vect3 grads = P1Vector(points[0], points[1], points[2], j);
                for(int i=0; i<3; i++)
                    A(3*t+i, trg[j]) = grads(i);
            }
        }
        return A;
    }

    /**
     * Computes the total solid angle of a surface for a point p and tells whether p is inside the mesh or not.
     **/
    /*
    bool Mesh::contains_point(const Vect3& p) const {

        double solangle = 0.0;
        for (int itrg=0;itrg<ntrgs;itrg++)
            solangle += p.solangl(pts[trgs[itrg][0]],pts[trgs[itrg][1]], pts[trgs[itrg][2]]);

        if (std::abs(solangle) < 1e3*std::numeric_limits<double>::epsilon()) {
            return false;
        } else if (std::abs(solangle + 4*M_PI) < 1e3*std::numeric_limits<double>::epsilon()) {
            return true;
        } else if (std::abs(solangle - 4*M_PI) < 1e3*std::numeric_limits<double>::epsilon()) {
            std::cerr << "Mesh::contains_point(" << p << ") Error. Are you sure the mesh is properly oriented?\n";
            return false;
        } else {
            std::cerr << "Mesh::contains_point(" << p << ") Error. Are you sure the mesh is closed?\n"
                      << std::abs(solangle) << std::endl;
            exit(1);
        }
    } */

    Vector Mesh::areas() const {
        Vector my_areas(nbTrgs());
        for(int i = 0; i < nbTrgs(); ++i)
            my_areas(i) = triangle(i).getArea();

        return my_areas;
    }

    void Mesh::update_triangles() {
        for(int i = 0; i < ntrgs; ++i) {
            trgs[i].normal() = pts[trgs[i][0]].normal( pts[trgs[i][1]] , pts[trgs[i][2]] ); // TODO: not Normalized
            trgs[i].area() = trgs[i].normal().norm() / 2.0;
        }
    }

    void Mesh::recompute_normals() {
        for (int p = 0; p < nbPts(); ++p) {
            Vect3 my_normal(0);
            for (intSet::iterator it = links[p].begin(); it != links[p].end(); ++it) {
                my_normal += trgs[*it].normal().normalize();
            }
            my_normal.normalize();
            normals[p] = my_normal;
        }
    }

    bool Mesh::has_self_intersection() const {
        bool selfIntersects = false;
        for (int i = 0; i < ntrgs; ++i) {
            const Triangle& T1 = triangle(i);
            for (int j = i+1; j < ntrgs; ++j) {
                const Triangle& T2 = triangle(j);
                if (!T1.contains(T2.s1()) && !T1.contains(T2.s2()) && !T1.contains(T2.s3())) {
                    if (triangle_intersection(*this, i, *this, j)) {
                        // selfIntersects = true;
                        std::cout << "triangles " << i << " and " << j << " are intersecting" <<std::endl;
                    }
                }
            }
        }
        return selfIntersects;
    }

    bool Mesh::intersection(const Mesh& m) const {
        bool intersects = false;
        for(int i = 0; i < ntrgs; ++i) {
            for(int j = 0; j < m.nbTrgs(); ++j) {
                intersects = intersects | triangle_intersection(*this, i, m, j);
            }
        }
        return intersects;
    }

    bool Mesh::has_correct_orientation() const {
        // First : check the local orientation (that all the triangles are all oriented in the same way)
        // define the triangle edges as (first point, second point)
        // if a triangle edge is ordered with (lower index, higher index) keep it in a edge_list
        // if not, exchange the two indices and put in a flipped_edge_list (lower index, higher index)
        // Transform the edge_list and flipped_edge_list unambigously into lists of numbers
        std::list<int> edge_list;
        std::list<int> flipped_edge_list;
        int radix = 10^(int(log2((double)npts)/log2(10.0))+1);
        for(int i=0;i<ntrgs;++i) {
            for(int j=1;j<=3;++j) {
                if (trgs[i].som(j)<trgs[i].next(j))
                    edge_list.push_back(trgs[i].som(j)*radix+trgs[i].next(j));
                else
                    flipped_edge_list.push_back(trgs[i].next(j)*radix+trgs[i].som(j));
            }
        }
        // Sort these two lists: they should be identical: if not, there is a local orientation problem.
        edge_list.sort();
        flipped_edge_list.sort();

        // check the global orientation: that the triangles are correctly oriented (outward-pointing normal)
        // compute the bounding box:
        double xmax = std::numeric_limits<int>::min();
        double ymax = std::numeric_limits<int>::min();
        double zmax = std::numeric_limits<int>::min();
        double xmin = std::numeric_limits<int>::max();
        double ymin = std::numeric_limits<int>::max();
        double zmin = std::numeric_limits<int>::max();
        for(int i=0; i<npts; ++i) {
            xmin = std::min(xmin, point(i).x());
            ymin = std::min(ymin, point(i).y());
            zmin = std::min(zmin, point(i).z());
            xmax = std::max(xmax, point(i).x());
            ymax = std::max(ymax, point(i).y());
            zmax = std::max(zmax, point(i).z());
        }
        Vect3 bbmin(xmin, ymin, zmin);
        Vect3 bbmax(xmax, ymax, zmax);
        Vect3 bbcenter = 0.5 * (bbmin + bbmax);

        // check the center of the bounding box is inside the mesh
        bool in_mesh = contains_point(bbcenter);
        if (!in_mesh)
            std::cerr << "Global orientation problem..." << std::endl << std::endl;

        if (flipped_edge_list != edge_list)
            std::cerr << "Local orientation problem..." << std::endl << std::endl;

        return in_mesh && (flipped_edge_list == edge_list);
    }

} // end namespace OpenMeeg
