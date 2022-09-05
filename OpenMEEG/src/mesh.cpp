// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <sstream>
#include <stack>
#include <algorithm>

#include <constants.h>
#include <mesh.h>
#include <MeshIO.h>
#include <geometry.h>

namespace OpenMEEG {

    //  We need a shared_ptr TODO

    Geometry* Mesh::create_geometry(Geometry* geom) {
        if (geom!=nullptr)
            return geom;
        return new Geometry;
    }

    Mesh::Mesh(const unsigned nv,const unsigned nt,Geometry* geometry): Mesh(geometry) {
        geom->vertices().reserve(nv);
        triangles().reserve(nt);
    }

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
        if (current_barrier())
            std::cout << "\t\tCurrent barrier : " << current_barrier() << std::endl;
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

    TriangleIndices Mesh::triangle(const Triangle& t) const {
        auto ind = [&](const unsigned i) { return &t.vertex(i)-&geometry().vertices()[0]; };
        return TriangleIndices(ind(0),ind(1),ind(2));
    }

    Triangle& Mesh::add_triangle(const TriangleIndices inds) {
        Vertex* vptrs[3];
        for (unsigned i=0; i<3; ++i)
            vptrs[i] = &geometry().vertices().at(inds[i]);
        triangles().push_back(vptrs);
        return triangles().back();
    }

    void Mesh::reference_vertices(const IndexMap& indmap) {
        for (const auto& mapping : indmap)
            vertices().push_back(&geometry().vertices().at(mapping.second));
    }

    void Mesh::clear() {
        vertices().clear();
        triangles().clear();
        mesh_name.clear();
        vertex_triangles.clear();
        outermost_ = false;
    }

    /// Update triangles area/normal, update vertex triangles and vertices normals if needed

    void Mesh::update(const bool topology_changed) {

        // If indices are not set, we generate them for sorting edge and testing orientation

        if (topology_changed) {
            make_adjacencies();
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
        Normal N(0);
        for (const auto& triangle : triangles(v))
            N += triangle->normal();
        N.normalize();
        return N;
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
        update(true);
    }

    /// Smooth Mesh

    void Mesh::smooth(const double& smoothing_intensity,const unsigned& niter) {

        typedef std::map<const Vertex*,std::set<Vertex>> Neighbors;
        Neighbors neighbors;
        for (const auto& vertex : vertices())
            for (const auto& triangle : triangles(*vertex))
                for (unsigned k=0;k<3;++k)
                    if (&triangle->vertex(k)!=vertex)
                        neighbors[vertex].insert(triangle->vertex(k));


        for (unsigned n=0; n<niter; ++n) {
            Vertices new_pts(vertices().size());
            for (const auto& vertex : vertices()) {
                Vertex newpt;
                const auto& vneighbors = neighbors.at(vertex);
                for (const auto& neighbor : vneighbors)
                    newpt += (smoothing_intensity*(neighbor-*vertex));
                new_pts.push_back(newpt/vneighbors.size());
            }
            unsigned i = 0;
            for (auto& vertex : vertices())
                *vertex = new_pts[i++];
        }
        update(false); // Updating triangles (areas + normals)
    }

    /// Sq. Norm Surface Gradient: square norm of the surfacic gradient of the P1 and P0 elements

    void Mesh::gradient_norm2(SymMatrix &A) const {

        /// Vertices (diagonal of the matrix).

        for (const auto& vp : vertices()) {
            const Vertex& v1 = *vp;
            const unsigned index = v1.index();
            for (const auto& t : triangles(v1)) {
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
                    if (tit1->intersects(*tit2)) {
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
                intersects = intersects | triangle1.intersects(triangle2);
        return intersects;
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
        update(true);
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

    void Mesh::check_consistency(const std::string& when) const {
        // check that all vertices lead to triangles whose edges are defined
        for (auto& V1 : vertices()) {
            for (auto& tp1 : triangles(*V1)) {
                try {
                    tp1->edge(*V1);
                } catch (const OpenMEEG::UnknownVertex&) {
                    std::ostringstream oss;
                    oss << "Mesh " << name() << " invalid    during " << when << ", requested triangle vertex address:" << std::endl << "  " << &V1 << " (" << V1 << ")" << std::endl << "but valid triangle vertex addresses are:" << std::endl;
                    for (unsigned i=0;i<3;++i)
                        oss << "  " << &(tp1->vertex(i)) << std::endl;
                    throw OpenMEEG::UnknownVertex(oss.str());
                }
            }
        }
        std::cout << "Mesh " << name() << " consistent during " << when << std::endl;
    }
}
