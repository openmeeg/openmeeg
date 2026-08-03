// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <cgal_lib.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>

#include <CGAL/Polyhedron_incremental_builder_3.h>

// Simplification function

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

namespace OpenMEEG {

    using SurfaceMesh = CGAL::Surface_mesh<K::Point_3>;
    using CGALDomain  = CGAL::Polyhedral_mesh_domain_3<SurfaceMesh,K>;

    // Sizing field

    struct Planes {

        Planes(const OpenMEEG::Matrix& M,const double field_size): mat(M),fs(field_size) { }

        FT operator()(const K::Point_3& p,const int,const CGALDomain::Index&) const {
            bool inside = true;
            for (unsigned i=0; i<mat.nlin(); ++i) {
                const OpenMEEG::Vect3 v(p.x()-mat(i,0),p.y()-mat(i,1),p.z()-mat(i,2));
                if (v(0)*mat(i,3)+v(1)*mat(i,4)+v(2)*mat(i,5)<0.0)
                    inside = false;
            }
            return (inside) ? fs/3.0 : fs;
        }

    private:

        const OpenMEEG::Matrix mat;
        double fs;
    };

    struct Spheres {

        Spheres(const OpenMEEG::Matrix& M,const double field_size): mat(M),fs(field_size) { }

        FT operator()(const K::Point_3& p,const int,const CGALDomain::Index&) const {
            bool inside = false;
            for (unsigned i=0; i<mat.nlin(); ++i) {
                const OpenMEEG::Vect3 v(p.x()-mat(i,0),p.y()-mat(i,1),p.z()-mat(i,2));
                if (v.norm()<mat(i,3))
                    inside = true;
            }
            return (inside) ? fs/3.0 : fs;
        }

    private:

        const OpenMEEG::Matrix mat;
        double fs;
    };

    SurfaceMesh Mesh2SurfaceMesh(const Mesh& m) {

        SurfaceMesh sm;

        std::map<const Vertex*,SurfaceMesh::Vertex_index> vmap;

        for (auto vit=m.vertices().begin(); vit!=m.vertices().end(); ++vit)
            vmap[*vit] = sm.add_vertex(K::Point_3((*vit)->x(),(*vit)->y(),(*vit)->z()));

        for (auto tit=m.triangles().begin(); tit!=m.triangles().end(); ++tit)
            sm.add_face(vmap[&tit->vertex(0)],vmap[&tit->vertex(1)],vmap[&tit->vertex(2)]);

        return sm;
    }

    /// Refine a mesh.

    Mesh cgal_refine(const Mesh& mesh,const double radius_bound,const double distance_bound,const std::string& sizing_field) {

        using Triangulation = CGAL::Mesh_triangulation_3<CGALDomain>::type;
        using C3t3P         = CGAL::Mesh_complex_3_in_triangulation_3<Triangulation>;
        using Criteria      = CGAL::Mesh_criteria_3<Triangulation>;

        // Create input mesh and domain.

        SurfaceMesh sm = Mesh2SurfaceMesh(mesh);
        CGALDomain  domain(sm);

        // Create the refined mesh.

        C3t3P c3t3;

        if (sizing_field.empty()) {

            Criteria criteria(CGAL::parameters::facet_angle(30).facet_size(radius_bound).facet_distance(distance_bound));
            c3t3 = CGAL::make_mesh_3<C3t3P>(domain,criteria,CGAL::parameters::no_exude(),CGAL::parameters::no_perturb());

        } else {

            Matrix field(sizing_field);
            if (field.ncol()==6) {
                Planes planes(field,radius_bound);
                Criteria criteria(CGAL::parameters::facet_angle(30).facet_size(planes).facet_distance(distance_bound));
                c3t3 = CGAL::make_mesh_3<C3t3P>(domain,criteria,CGAL::parameters::no_exude(),CGAL::parameters::no_perturb());
            } else if (field.ncol()==4) {
                Spheres spheres(field,radius_bound);
                Criteria criteria(CGAL::parameters::facet_angle(30).facet_size(spheres).facet_distance(distance_bound));
                c3t3 = CGAL::make_mesh_3<C3t3P>(domain,criteria,CGAL::parameters::no_exude(),CGAL::parameters::no_perturb());
            } else {
                std::cerr << "Error: file should contain either 4 or 6 columns" << std::endl;
            }
        }

        return CGAL_to_OM(c3t3);
    }

    /// Decimate safely a mesh.

    Mesh cgal_decimate(const Mesh& mesh,const unsigned nb_points) {

        namespace SMS = CGAL::Surface_mesh_simplification;

        // Create input mesh

        SurfaceMesh sm = Mesh2SurfaceMesh(mesh);

        const unsigned nb_edges = 3*nb_points-6;

        std::cout << "Target: less than " << nb_points << " vertices <=> " << nb_edges << " edges.\n";

        SMS::Edge_count_stop_predicate<SurfaceMesh> stop(nb_edges);

        const int nb_edge_removed = SMS::edge_collapse(sm,stop,CGAL::parameters::get_cost(SMS::Edge_length_cost<SurfaceMesh>()).get_placement(SMS::Midpoint_placement<SurfaceMesh>()));

        std::cout << "Number of vertices  = " << sm.number_of_vertices() << std::endl
                  << "Number of triangles = " << sm.number_of_faces()    << std::endl;

        // write the output mesh

        Mesh m(sm.number_of_vertices(),sm.number_of_faces());
        Geometry& geom = m.geometry();

        // Insert points in the geometry.

        std::map<SurfaceMesh::Vertex_index,unsigned> vertex_index_mapping;
        for (SurfaceMesh::Vertex_index v : sm.vertices()) {
            const K::Point_3& p     = sm.point(v);
            const unsigned    index = geom.add_vertex(Vertex(CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z())));
            vertex_index_mapping[v] = index;
            m.vertices().push_back(&geom.vertices()[index]);
        }

        for (SurfaceMesh::Face_index f : sm.faces()) {

            unsigned k = 0;
            unsigned inds[3];
            for (auto v : CGAL::vertices_around_face(sm.halfedge(f),sm))
                inds[k++] = vertex_index_mapping[v];

            m.add_triangle(TriangleIndices(inds));
        }

        #if 0
        // Insert points to the mesh.

        std::set<const Vertex*> mesh_vertices;
        for (auto& triangle : m.triangles())
            for (auto& vertex : triangle)
                if (mesh_vertices.insert(vertex).second)
                    m.vertices().push_back(vertex);
        #endif

        m.update();
        m.correct_global_orientation();

        std::cout << "Finished...\n" << nb_edge_removed << " edges removed.\n" << m.vertices().size() << " final nb_points.\n";

        return m;
    }
}
