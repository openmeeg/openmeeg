// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

// for verbosity
#define CGAL_MESH_3_VERBOSE

#include <geometry.h>
#include <mesh.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/number_utils.h>

using K  = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = K::FT;

namespace OpenMEEG {

    //  Implicit Sphere and HemiSphere functions.

    class SphereFunction {
    public:

        SphereFunction(const FT sqradius): squared_radius(sqradius) { }

        FT operator()(const K::Point_3& p) const {
            return CGAL::square(p.x())+CGAL::square(p.y())+CGAL::square(p.z())-squared_radius;
        }

    private:

        FT squared_radius;
    };

    class HemiSphereFunction {
    public:

        HemiSphereFunction(const FT sqradius): squared_radius(sqradius) { }

        FT operator()(const K::Point_3& p) const {
            const double x2 = CGAL::square(p.x());
            const double y2 = CGAL::square(p.y());
            const double z2 = CGAL::square(p.z());
            const double dist_sphere = x2+y2+z2-squared_radius;
            if (p.z()>0) {
                return (dist_sphere>0) ? dist_sphere : std::max(dist_sphere,-z2);
            } else {
                return std::min(x2,squared_radius)+std::min(y2,squared_radius)+z2;
            }
        }

    private:

        FT squared_radius;
    };

    template <typename C3t3>
    Mesh CGAL_to_OM(const C3t3& c3t3) {

        Mesh m(c3t3.triangulation().number_of_vertices(),c3t3.number_of_facets());
        Geometry& geom = m.geometry();

        using Triangulation = C3t3::Triangulation;

        std::map<typename C3t3::Vertex_handle,unsigned> vertex_index_mapping;
        for (typename Triangulation::Finite_vertices_iterator vit = c3t3.triangulation().finite_vertices_begin(); vit!=c3t3.triangulation().finite_vertices_end(); ++vit) {
            const typename Triangulation::Point& p = vit->point();
            const unsigned index = geom.add_vertex(Vertex(CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z())));
            vertex_index_mapping[vit] = index;
            m.vertices().push_back(&geom.vertices()[index]);
        }

        for (typename C3t3::Facets_in_complex_iterator fit = c3t3.facets_in_complex_begin(); fit!=c3t3.facets_in_complex_end(); ++fit) {
            const typename Triangulation::Cell_handle cell = fit->first;
            const int& index = fit->second;
            const unsigned index1 = vertex_index_mapping[cell->vertex(c3t3.triangulation().vertex_triple_index(index,0))];
            const unsigned index2 = vertex_index_mapping[cell->vertex(c3t3.triangulation().vertex_triple_index(index,1))];
            const unsigned index3 = vertex_index_mapping[cell->vertex(c3t3.triangulation().vertex_triple_index(index,2))];
            m.add_triangle({ index1, index2, index3 });
        }

        m.update();
        m.info();

        m.correct_global_orientation();

        return m;
    }

    /// Create a mesh from an implicit function.

    template <typename Function,typename BoundingSphere>
    Mesh cgal_mesh_function(const Function& f,const BoundingSphere& bound,const double radius_bound,const double distance_bound) {

        using MeshDomain    = CGAL::Labeled_mesh_domain_3<K>;
        using Triangulation = CGAL::Mesh_triangulation_3<MeshDomain,CGAL::Default,CGAL::Parallel_if_available_tag>::type;
        using C3t3S         = CGAL::Mesh_complex_3_in_triangulation_3<Triangulation>;
        using Criteria      = CGAL::Mesh_criteria_3<Triangulation>;

        // Defining the domain,

        const MeshDomain& domain = MeshDomain::create_implicit_mesh_domain(f,bound,CGAL::parameters::relative_error_bound(1e-6));

        // Mesh criteria

        Criteria crit(CGAL::parameters::facet_angle(30).facet_size(radius_bound).facet_distance(distance_bound));

        // meshing domain

        const C3t3S& c3t3 = CGAL::make_mesh_3<C3t3S>(domain,crit,CGAL::parameters::no_exude(),CGAL::parameters::no_perturb());

# if 0
        // If you want want to add initial points on the hemisphere circle (for a better definition),
        // have a look here (it probably needs to construct the facets also ).

        std::pair<Tr::Point,unsigned> p[init_points];
        for ( unsigned iip = 0; iip < init_points; ++iip) {
            p[iip] = std::make_pair(Tr::Point(hemisphere_radius*std::cos(2.*Pi/init_points*iip),hemisphere_radius*std::sin(2.*Pi/init_points*iip),0),0);
        }
        c3t3.insert_surface_points(&p[0],&p[init_points-1]);
        CGAL::refine_mesh_3<C3t3>(c3t3,hdomain,criteria,no_exude(),no_perturb());
#endif

        return CGAL_to_OM(c3t3);
    }

    Mesh cgal_refine(const Mesh& m_in,const double radius_bound,const double distance_bound,const std::string& sizing_field);
    Mesh cgal_decimate(const Mesh& m_in,const unsigned nb_points);

    #ifdef CGAL_ImageIO
    Mesh cgal_mesh_3Dlabeled_image(const char* input_filename,const double radius_bound,const double distance_bound);
    Mesh cgal_mesh_3Dlevelset_image(const char* input_filename,const double levelset_value,const bool positive_inside,const double radius_bound,const double distance_bound);
    #endif
}
