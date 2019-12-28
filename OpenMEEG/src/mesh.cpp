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

#include <sstream>
#include <stack>
#include <algorithm>

#include <mesh.h>
#include <mesh_ios.h>
#include <geometry.h>

#include <Triangle_triangle_intersection.h>

namespace OpenMEEG {

    Mesh::Mesh(): geom(new Geometry),mesh_triangles() { }

    Mesh::Mesh(const unsigned nv,const unsigned nt): geom(new Geometry) {
        geom->vertices().reserve(nv);
        triangles().reserve(nt);
    }

    Mesh::Mesh(Vertices&& vs,Triangles&& ts): geom(new Geometry) {
        geom->vertices() = vs;
        triangles()      = ts;
    }

    Mesh::Mesh(const std::string& filename,const bool verbose,const std::string& name):
        mesh_name(name),geom(new Geometry)
    {
        load(filename,verbose);
    }

    #if 0
    Mesh::Mesh(Geometry& g,const std::string name): mesh_name(name),geom(&g) { }

    void Mesh::copy(const Mesh& m) {

        if (m.all_vertices_) {
            all_vertices_ = std::make_shared<Vertices>();
            all_vertices_->reserve(m.vertices().size());
            std::map<const Vertex*,Vertex*> map; // for the triangles
            for (const auto& vertex : m.vertices()) {
                add_vertex(*vertex);
                map[vertex] = vertices().back();
            }
            for (const auto& triangle : m.triangles()) {
                Triangle t(map[&triangle.vertex(0)],map[&triangle.vertex(1)],map[&triangle.vertex(2)]);
                triangles().push_back(t);
            }
            update();
        } else {
            all_vertices_ = m.all_vertices_;
            for (const auto& triangle : m.triangles())
                triangles().push_back(triangle);
            build_mesh_vertices();
        }
        outermost_ = m.outermost_;
        mesh_name  = m.mesh_name;
    }
    #endif

    /// Print informations about the mesh 

    void Mesh::info(const bool verbose) const {
        std::cout << "Info:: Mesh name/ID : "  << name() << std::endl;
        std::cout << "\t\t# vertices  : " << vertices().size() << std::endl;
        std::cout << "\t\t# triangles : " << triangles().size() << std::endl;
        std::cout << "\t\tEuler characteristic : " << vertices().size()-3*triangles().size()/2+triangles().size() << std::endl;

        double min_area = std::numeric_limits<double>::max();
        double max_area = 0.;
        for (const auto& triangle : triangles()) {
            min_area = std::min(triangle.area(),min_area);
            max_area = std::max(triangle.area(),max_area);
        }
        std::cout << "\t\tMin Area : " << min_area << std::endl;
        std::cout << "\t\tMax Area : " << max_area << std::endl;
        if (verbose) {
            std::cout << "Indices :" << std::endl;
            for (const auto& vertex : vertices())
                std::cout << "[" << *vertex << "] = " << vertex->index() << std::endl;
            for (const auto& triangle : triangles())
                std::cout << "[[" << triangle.vertex(0)
                          << "] , [" << triangle.vertex(1)
                          << "] , ["<< triangle.vertex(2)
                          << "]] \t = " << triangle.index() << std::endl;
        }
    }

    void Mesh::build_mesh_vertices() {

        // Sets do not preserve the order, and we would like to preserve it so we push_back in the vector as soon as the element is unique.

        std::set<const Vertex *> mesh_v;
        vertices().clear();
        for (auto& triangle : triangles())
            for (auto& vertex : triangle)
                if (mesh_v.insert(vertex).second)
                    vertices().push_back(vertex);
    }

    void Mesh::clear() {
        std::cerr << geom << std::endl;
        std::cerr << "1" << std::endl;
        vertices().clear();
        triangles().clear();
        std::cerr << "3" << std::endl;
        mesh_name.clear();
        std::cerr << "6" << std::endl;
        links_.clear();
        std::cerr << "7" << std::endl;
        outermost_ = false;
        std::cerr << "8" << std::endl;
        std::cerr << "Destroyed" << std::endl;
    }

#if 0
    /// Add vertex to the mesh/geometry.

    void Mesh::add_vertex(const Vertex& V) {

        geom->add_vertex(V);

            if (std::find(vertices().begin(),vertices().end(),&(*vit))==vertices().end())
                vertices().push_back(&(*vit));
        }
    }
#endif

    /// Update triangles area/normal, update links and vertices normals if needed
    void Mesh::update() {

        make_adjacencies();

        // If indices are not set, we generate them for sorting edge and testing orientation

        if (!finalized) {
            finalized = true;
            if (vertices().front()->index()==unsigned(-1))
                generate_indices();
            correct_local_orientation();
        }

        // Compute triangles' normals and areas (after having the mesh locally reoriented)

        for (auto& triangle : triangles()) {
            Vect3 normaldir   = crossprod(triangle.vertex(0)-triangle.vertex(1),triangle.vertex(0)-triangle.vertex(2));
            triangle.area()   = normaldir.norm()/2.0;
            triangle.normal() = normaldir.normalize();
        }
    }

    /// Compute normals at vertices.

    Normal Mesh::normal(const Vertex& v) const {
        Normal _normal(0);
        for (TrianglesRefs::const_iterator tit = links_.at(&v).begin(); tit != links_.at(&v).end(); ++tit)
            _normal += (*tit)->normal();
        _normal.normalize();
        return _normal;
    }

    /// Add a mesh (assumes geometry points are of sufficient size.

    void Mesh::add_mesh(const Mesh& m) {
        std::map<const Vertex*,Vertex*> vmap;
        for (const auto& vertex : m.vertices())
            vmap[vertex] = &geom->vertices().at(geom->add_vertex(*vertex));
        auto vertex = [&](const Triangle& t,const unsigned ind) { return vmap.at(&t.vertex(ind)); };
        for (const auto& triangle : m.triangles())
            triangles().push_back(Triangle(vertex(triangle,0),vertex(triangle,1),vertex(triangle,2)));
    }

    /// Merge two meshes into one (without duplicating vertices).

    void Mesh::merge(const Mesh& m1, const Mesh& m2) {
        clear();
        geom->vertices().reserve(m1.vertices().size()+m2.vertices().size());
        add_mesh(m1);
        add_mesh(m2);
        update();
    }

    /// Smooth Mesh

    void Mesh::smooth(const double& smoothing_intensity,const unsigned& niter) {

        typedef std::vector<std::set<Vertex>> Neighbors;
        Neighbors neighbors(vertices().size());
        unsigned i = 0;
        for (const auto& vertex : vertices()) {
            for (const auto& triangle : adjacent_triangles(*vertex))
                for (unsigned  k=0;k<3;++k)
                    if (triangle->vertex(k)==*vertex) // TODO: Is a pointer comparison sufficient ?
                        neighbors[i].insert(triangle->vertex(k));
            ++i; // Can't we use the vertex index instead ? TODO
        }

        Vertices new_pts(vertices().size());
        for (unsigned n = 0; n < niter; ++n) {
            Neighbors::const_iterator nit = neighbors.begin();
            for (const auto& vertex : vertices()) {
                new_pts.push_back(*vertex);
                auto& newpt = new_pts.back();
                for (std::set<Vertex>::const_iterator it = nit->begin();it!=nit->end();++it)
                    newpt += (smoothing_intensity*(*it-*vertex))/nit->size();
                ++nit;
            }
            for (auto& vertex : vertices())
                *vertex = new_pts[i];
            new_pts.clear();
        }
        update(); // Updating triangles (areas + normals)
    }

    /// Sq. Norm Surface Gradient: square norm of the surfacic gradient of the P1 and P0 elements

    void Mesh::gradient_norm2(SymMatrix &A) const {

        /// Vertices (diagonal of the matrix).

        for (const auto& vp : vertices()) {
            const Vertex& v1 = *vp;
            const unsigned index = v1.index();
            for (const auto& t : links_.at(vp)) {
                const Edge& edge = t->edge(v1);
                const Vertex& v2 = edge.vertex(0);
                const Vertex& v3 = edge.vertex(1);
                A(index,index) += P1gradient(v1,v2,v3).norm2()*sqr(t->area());
            }
        }

        // Edges

        for (const auto& triangle: triangles()) {
            const double area2 = sqr(triangle.area());
            const Edges& edges = triangle.edges();
            for (unsigned i=0;i<3;++i) {
                const Vertex& V1 = edges[i].vertex(0);
                const Vertex& V2 = edges[i].vertex(1);

                const unsigned ind1 = V1.index();
                const unsigned ind2 = V2.index();

                //  This is a symmetric matrix. We only need to fill the lower half matrix.

                if (ind1<ind2) {
                    const Vertex& V3 = triangle.vertex(i);
                    A(ind1,ind2) += dotprod(P1gradient(V1,V2,V3),P1gradient(V2,V3,V1))*area2;
                }
            }
        }

        // P0 gradients: loop on triangles
        if (!outermost_) // if it is an outermost mesh: p=0 thus no need for computing it
            for (const auto& triangle1 : triangles()) {
                A(triangle1.index(),triangle1.index()) = 0.;
                for (const auto& triangle2 : adjacent_triangles(triangle1))
                    if (triangle1.index()<triangle2->index()) // sym matrix only lower half
                        A(triangle1.index(), triangle2->index()) += P0gradient_norm2(triangle1,*triangle2)*triangle1.area()*triangle2->area();
            }
    }

    /// Laplacian Mesh: good approximation of Laplace-Beltrami operator
    // "Discrete Laplace Operator on Meshed Surfaces". by Belkin, Sun, Wang

    void Mesh::laplacian(SymMatrix& A) const {
        for (const auto& triangle : triangles())
            for (const auto& edge : triangle.edges()) {
                const Vertex& V1 = edge.vertex(0);
                const Vertex& V2 = edge.vertex(1);
                const unsigned ind1 = V1.index();
                const unsigned ind2 = V2.index();
                if (ind1<ind2) { // sym matrix only lower half
                    const double h = (V2-V1).norm();
                    A(ind1,ind2) += -triangle.area()/(12*Pi*sqr(h))*exp(-sqr(h)/(4*h));
                }
            }

        for (const auto& vertex : vertices())
            A(vertex->index(),vertex->index()) = -A.getlin(vertex->index()).sum();
    }

    bool Mesh::has_self_intersection() const {

        bool selfIntersects = false;
        for (auto tit1=triangles().begin();tit1!=triangles().end();++tit1)
            for (auto tit2=tit1;tit2!=triangles().end();++tit2)
                if (!tit1->contains(tit2->vertex(0)) && !tit1->contains(tit2->vertex(1)) && !tit1->contains(tit1->vertex(2)))
                    if (triangle_intersection(*tit1, *tit2)) {
                        selfIntersects = true;
                        std::cout << "Triangles " << tit1->index() << " and " << tit2->index() << " are intersecting." << std::endl;
                    }
        return selfIntersects;
    }

    double Mesh::solid_angle(const Vect3& p) const {
        double solangle = 0.0;
        for (const auto& triangle : triangles())
            solangle += p.solid_angle(triangle.vertex(0),triangle.vertex(1),triangle.vertex(2));
        return solangle;
    }

    bool Mesh::intersection(const Mesh& m) const {
        bool intersects = false;
        for (const auto& triangle1 : triangles())
            for (const auto& triangle2 : m.triangles())
                intersects = intersects | triangle_intersection(triangle1,triangle2);
        return intersects;
    }

    bool Mesh::triangle_intersection(const Triangle& T1, const Triangle& T2) const {
        const Vect3& p1 = T1.vertex(0);
        const Vect3& q1 = T1.vertex(1);
        const Vect3& r1 = T1.vertex(2);
        const Vect3& p2 = T2.vertex(0);
        const Vect3& q2 = T2.vertex(1);
        const Vect3& r2 = T2.vertex(2);

        double pp1[3] = { p1.x(), p1.y(), p1.z() };
        double qq1[3] = { q1.x(), q1.y(), q1.z() };
        double rr1[3] = { r1.x(), r1.y(), r1.z() };
        double pp2[3] = { p2.x(), p2.y(), p2.z() };
        double qq2[3] = { q2.x(), q2.y(), q2.z() };
        double rr2[3] = { r2.x(), r2.y(), r2.z() };
        return tri_tri_overlap_test_3d(pp1,qq1,rr1,pp2,qq2,rr2);
    }

    const TrianglesRefs& Mesh::triangles(const Vertex& V) const {
        AdjacencyMap::const_iterator it = links_.find(&V);
        if (it!=links_.end())
            return it->second;

        // TODO ??? Throw an exception ?
        static TrianglesRefs a;
        return a;
    }

    void Mesh::load(const std::string& filename,const bool verbose) {
        clear();
        MeshIO* io = MeshIO::create(filename);
        if (verbose)
            std::cout << "loading : " << filename << " as a \"" << io->name() << "\" file."<< std::endl;
        io->load(*this);
        delete io;
        if (verbose)
            info();
        update();
        generate_indices();
    }

    void Mesh::generate_indices() {
        unsigned index = 0;
        for (auto& vertex : vertices())
            vertex->index() = index++;
        for (auto& triangle : triangles())
            triangle.index() = index++;
    }

    void Mesh::save(const std::string& filename) const {
        MeshIO* io = MeshIO::create(filename);
        io->save(*this);
        delete io;
    }

    const Mesh::EdgeMap Mesh::compute_edge_map() const {

        // Associate an integer with each edge.
        // Well oriented inner edges will be mapped to 0, border edges to 1 and badly oriented edges to 2.
        // The algorithm goes through each triangle edge e=(first vertex, second vertex)
        // If e is ordered with (lower index, higher index) add 1 to its map else remove 1.

        EdgeMap edgemap;
        auto lambda = [&](const Vertex& V1,const Vertex& V2,const int incr) {
            const auto& index = std::make_pair(&V1,&V2);
            if (edgemap.count(index)==0)
                edgemap[index] = incr;
            else
                edgemap[index] += incr;
        };

        for (const auto& triangle : triangles())
            for (const auto& edge : triangle.edges()) {
                const Vertex& V1 = edge.vertex(0);
                const Vertex& V2 = edge.vertex(1);
                if (V1.index()>V2.index())
                    lambda(V1,V2,1);
                else
                    lambda(V2,V1,-1);
            }

        return edgemap;
    }

    void Mesh::correct_local_orientation() {
        if (!has_correct_orientation()) {
            std::cerr << "Reorienting..." << std::endl << std::endl;
            std::stack<Triangle*>    triangle_stack;
            std::map<Triangle*,bool> reoriented_triangles;
            triangle_stack.push(&triangles().front());
            reoriented_triangles[&triangles().front()] = true;

            const auto has_same_edge = [](const Edges& edges1,const Edges& edges2) {
                for (const auto& edge2 : edges2)
                    for (const auto& edge1 : edges1)
                        if (edge1==edge2)
                            return true;
                return false;
            };

            while (!triangle_stack.empty()) {
                const Triangle& t1     = *triangle_stack.top();
                const Edges&    edges1 = t1.edges();
                triangle_stack.pop();
                for (const auto& tp : adjacent_triangles(t1))
                    if (reoriented_triangles.count(tp)==0) {
                        triangle_stack.push(tp);
                        Triangle& t2 = *tp;
                        const Edges& edges2 = t2.edges();
                        if (has_same_edge(edges1,edges2))
                            t2.change_orientation();
                        reoriented_triangles[tp] = true;
                    }
            }
        }
    }

    /// warning: a mesh may not be closed (as opposite to an interface)
    /// this function is mainly needed when meshing closed mesh with CGAL.

    void Mesh::correct_global_orientation() {

        Vect3 center(0.,0.,0.);
        for (const auto& vertex : vertices())
            center += *vertex;

        center /= vertices().size();
        const double solangle = solid_angle(center);
        if (!almost_equal(solangle,-4*Pi)) {
            if (almost_equal(solangle,0.0)) {
                std::cout << "Center point :" << center << " is on the mesh." << std::endl;
            } else if (almost_equal(solangle,4*Pi)) {
                change_orientation();
            } else {
                std::cout << "Not a closed mesh." << std::endl;
            }
        }
    }

    bool Mesh::has_correct_orientation() const {

        /// Check the local orientation (that all the triangles are all oriented in the same way)

        const EdgeMap mape = compute_edge_map();

        for (EdgeMap::const_iterator eit = mape.begin(); eit != mape.end(); ++eit)
            if (std::abs(eit->second) == 2) {
                std::cerr << "Local orientation problem..." << std::endl << std::endl;
                return false;
            }

        return true;
    }
}
