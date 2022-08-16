// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

// for verbosity
#define CGAL_MESH_3_VERBOSE

#include <mesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>
#include <CGAL/Image_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Implicit_mesh_domain_3.h>

#include <CGAL/config.h>
#include <CGAL/version.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::FT FT;
typedef CGAL::Polyhedron_3<K>  Polyhedron;
typedef Polyhedron::HalfedgeDS HDS;

namespace OpenMEEG {
    template <typename C3t3>
    Mesh CGAL_to_OM(C3t3 c3t3) {
        typedef typename C3t3::Triangulation Tr;
        Mesh m(c3t3.triangulation().number_of_vertices(), c3t3.number_of_facets());

        std::map<typename C3t3::Vertex_handle, unsigned> V;
        unsigned inum = 0;
        for ( typename Tr::Finite_vertices_iterator vit = c3t3.triangulation().finite_vertices_begin(), end = c3t3.triangulation().finite_vertices_end(); vit != end; ++vit)
        {
            const typename Tr::Point& p = vit->point();
            m.add_vertex(Vertex(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z())));
            V[vit] = inum++;
        }

        for ( typename C3t3::Facets_in_complex_iterator fit = c3t3.facets_in_complex_begin(); fit != c3t3.facets_in_complex_end(); ++fit) {
            const typename Tr::Cell_handle cell = fit->first;
            const int& index = fit->second;
            const int index1 = V[cell->vertex(c3t3.triangulation().vertex_triple_index(index, 0))];
            const int index2 = V[cell->vertex(c3t3.triangulation().vertex_triple_index(index, 1))];
            const int index3 = V[cell->vertex(c3t3.triangulation().vertex_triple_index(index, 2))];
            Triangle t(m.vertices()[index1], m.vertices()[index2], m.vertices()[index3]);
            m.push_back(t);
        }

        m.update();
        m.correct_global_orientation();
        return m;
    }


    Mesh cgal_mesh_function(double sphere_radius, double hemisphere, double radius_bound, double distance_bound);

    #if CGAL_VERSION_NR >= CGAL_VERSION_NUMBER(4,6,0)
    Mesh cgal_refine(const Mesh& m_in, double radius_bound, double distance_bound, const char* sizing_field);

    Mesh cgal_decimate(const Mesh& m_in, unsigned nb_points);
    #endif

    #ifdef CGAL_ImageIO
    Mesh cgal_mesh_3Dlabeled_image(const char* input_filename, double radius_bound, double distance_bound);
    Mesh cgal_mesh_3Dlevelset_image(const char* input_filename,  double levelset_value, bool positive_inside, double radius_bound, double distance_bound);
    #endif
}
