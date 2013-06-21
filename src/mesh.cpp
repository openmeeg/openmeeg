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

#include "mesh.h"
#include "Triangle_triangle_intersection.h"

namespace OpenMEEG {

    // Mesh::Mesh(const Mesh& M) { *this = M; }


    // Mesh& Mesh::operator= (const Mesh& M) {
        // if (this != &M)
            // copy(M);
        // return *this;
    // }

    // void Mesh::copy(const Mesh& M) {
        // destroy();
        // if (M.npts != 0) {
            // _all_vertices = M.all_vertices();
            // ntrgs   = M.ntrgs;
            // pts     = new Vect3[npts];
            // trgs    = new Triangle[ntrgs];
            // links   = new intSet[npts];
            // normals = new Vect3[npts];
            // for (int i=0; i < npts; i++) {
                // pts[i]     = M.pts[i];
                // links[i]   = M.links[i];
                // normals[i] = M.normals[i];
            // }
            // for (int i=0; i < ntrgs; i++)
                // trgs[i] = M.trgs[i];
        // }
    // }

    // void Mesh::destroy() {
        // if (npts != 0) {
            // delete [] pts;
            // delete [] trgs;
            // delete [] links;
            // delete [] normals;

            // npts = 0;
            // ntrgs = 0;
        // }
    // }

    void Mesh::update() {
        SetPVertex vset;
        for (Triangles::iterator tit = this->begin(); tit != this->end(); tit++) {
            tit->area()   = tit->normal().norm() / 2.0;
            vset.insert(&tit->s1().vertex());
            vset.insert(&tit->s2().vertex());
            vset.insert(&tit->s3().vertex());
        }
        vertices().insert(vertex_begin(), vset.begin(), vset.end()); // copy the set of mesh vertices into a vector
        vset.clear();
        links().reserve(vertices().size());
        size_t i = 0;
        for (const_vertex_iterator vit = vertex_begin(); vit != vertex_end(); vit++, i++) {
            for (const_iterator tit = this->begin(); tit != this->end(); tit++) {
                if (tit->contains(**vit)) {
                    links()[i].insert(*tit); 
                }
            }
        }
    }

    /**
     * Print informations about the mesh
    **/
    void Mesh::info() const {

        std::cout << "Mesh Info : "     << std::endl;
        std::cout << "\tName : "        << name() << std::endl;
        std::cout << "\t# points : "    << nb_vertices() << std::endl;
        std::cout << "\t# triangles : " << nb_triangles() << std::endl;
        std::cout << "\tEuler characteristic : " << nb_vertices() - 3.*nb_triangles()/2. + nb_triangles() << std::endl;

        double min_area = 0.;
        double max_area = 0.;
        for (const_iterator tit = this->begin(); tit != this->end(); tit++) {
            min_area = std::min(tit->area(), min_area);
            max_area = std::max(tit->area(), max_area);
        }
        std::cout << "\tMin Area : " << min_area << std::endl;
        std::cout << "\tMax Area : " << max_area << std::endl;
    }

    // void Mesh::update_triangles() {
        // for (const_iterator tit = this->begin(); tit != this->end(); tit++) {
            // tit->normal() = (*tit)(0)->vertex().normal( (*tit)(1) , (*tit)(2) ); // TODO: not Normalized
            // tit->area()   = tit->normal().norm() / 2.0;
        // }
    // }

    bool Mesh::has_self_intersection() const {
        bool selfIntersects = false;
        for (const_iterator tit1 = this->begin(); tit1 != this->end(); tit1++) {
            for (const_iterator tit2 = tit1; tit2 != this->end(); tit2++) {
                if (!tit1->contains(tit2->s1().vertex()) && !tit1->contains(tit2->s2().vertex()) && !tit1->contains(tit1->s3().vertex())) {
                    if (triangle_intersection(*tit1, *tit2)) {
                        selfIntersects = true;
                        std::cout << "triangles " << tit1->index() << " and " << tit2->index() << " are intersecting" << std::endl;
                    }
                }
            }
        }
        return selfIntersects;
    }

    bool Mesh::intersection(const Mesh& m) const {
        bool intersects = false;
        for (const_iterator tit1 = this->begin(); tit1 != this->end(); tit1++) {
            for (const_iterator tit2 = m.begin(); tit2 != m.end(); tit2++) {
                intersects = intersects | triangle_intersection(*tit1, *tit2);
            }
        }
        return intersects;
    }

    bool Mesh::triangle_intersection(const Triangle& T1, const Triangle& T2 ) const {
        const Vect3& p1 = T1.s1().vertex();
        const Vect3& q1 = T1.s2().vertex();
        const Vect3& r1 = T1.s3().vertex();
        const Vect3& p2 = T2.s1().vertex();
        const Vect3& q2 = T2.s2().vertex();
        const Vect3& r2 = T2.s3().vertex();

        double pp1[3] = {p1.x(), p1.y(), p1.z()};
        double qq1[3] = {q1.x(), q1.y(), q1.z()};
        double rr1[3] = {r1.x(), r1.y(), r1.z()};
        double pp2[3] = {p2.x(), p2.y(), p2.z()};
        double qq2[3] = {q2.x(), q2.y(), q2.z()};
        double rr2[3] = {r2.x(), r2.y(), r2.z()};
        return tri_tri_overlap_test_3d(pp1, qq1, rr1, pp2, qq2, rr2);
    }

    /* TODO
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
    } */

} // end namespace OpenMeeg

