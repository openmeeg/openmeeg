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

#ifndef OPENMEEG_CGAL_MESH_INCLUDE_H
#define OPENMEEG_CGAL_MESH_INCLUDE_H

// for verbosity
#define CGAL_MESH_3_VERBOSE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Gray_level_image_3.h> // TODO
#include <CGAL/Labeled_image_mesh_domain_3.h> // TODO
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_3::Robust_intersection_traits_3<K> Geom_traits;
// typedef CGAL::Polyhedron_3<K> Polyhedron;
// Domain 
typedef CGAL::Gray_level_image_3<Geom_traits::FT, Geom_traits::Point_3> Gray_level_image; // TODO
typedef CGAL::Implicit_surface_3<Geom_traits, Gray_level_image> GreySurface_3; // TODO
typedef CGAL::Labeled_image_mesh_domain_3<CGAL::Image_3,K> Mesh_domain; // TODO

typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

typedef Tr::Vertex_handle Vertex_handle;
typedef Tr::Geom_traits GT;


#if 0
typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Polyhedral_mesh_domain;

// Triangulation
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Vertex_handle Vertex_handle;

// from polygons
namespace CGAL {
    template <class Polyhedron, class Kernel, class Dummy_kernel> class AABB_polyhedral_oracle;
}

// sphere function
// for spheres
template <typename FT, typename Point>
class FT_to_point_Sphere_wrapper : public std::unary_function<Point, FT>
{
    double sqrd;
    public:
    FT_to_point_Sphere_wrapper(FT sqrd_) : sqrd(sqrd_) {}
    FT operator()(Point p) const { return (std::pow(p.x(),2)+std::pow(p.y(),2)+std::pow(p.z(),2)-sqrd); }
};
typedef FT_to_point_Sphere_wrapper<K::FT, K::Point_3> SphereFunction;
typedef CGAL::Implicit_surface_3<Geom_traits, SphereFunction> SphereSurface_3;

// Hemisphere function
template <typename FT, typename Point>
class FT_to_point_HemiSphere_wrapper : public std::unary_function<Point, FT>
{
    double sqrd;
    public:
    FT_to_point_HemiSphere_wrapper(FT sqrd_) : sqrd(sqrd_) {}
    FT operator()(Point p) const {

        double d_sphere = (std::pow(p.x(), 2) + std::pow(p.y(), 2) + std::pow(p.z(), 2) - sqrd); 
        double d_total;

        if ( p.z() > 0 ) {
            d_total = ( d_sphere > 0 ) ? d_sphere : -std::min(-d_sphere, std::pow(p.z(), 2));
        } else {
            d_total = std::min(std::pow(p.x(), 2), sqrd ) + std::min(std::pow(p.y(), 2), sqrd ) + std::pow(p.z(), 2);
        }

        return d_total;
    }
};
typedef FT_to_point_HemiSphere_wrapper<K::FT, K::Point_3> HemiSphereFunction;
typedef CGAL::Implicit_surface_3<Geom_traits, HemiSphereFunction> HemiSphereSurface_3;
#endif

OpenMEEG::Mesh CGAL_to_OM(C3t3 c3t3) {
    OpenMEEG::Mesh m(c3t3.triangulation().number_of_vertices(), c3t3.number_of_facets());

    std::map<Vertex_handle, unsigned> V;
    unsigned inum = 0;
    for ( Tr::Finite_vertices_iterator vit = c3t3.triangulation().finite_vertices_begin(), end = c3t3.triangulation().finite_vertices_end(); vit != end; ++vit)
    {
        const Point_3& p = vit->point();
        m.add_vertex(OpenMEEG::Vertex(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z())));
        V[vit] = inum++;
    }

    for ( Tr::Finite_facets_iterator fit = c3t3.triangulation().finite_facets_begin(); fit != c3t3.triangulation().finite_facets_end(); ++fit) {
        const Tr::Cell_handle cell = fit->first;
        const int& index = fit->second;
        const int index1 = V[cell->vertex(c3t3.triangulation().vertex_triple_index(index, 0))];
        const int index2 = V[cell->vertex(c3t3.triangulation().vertex_triple_index(index, 1))];
        const int index3 = V[cell->vertex(c3t3.triangulation().vertex_triple_index(index, 2))];
        OpenMEEG::Triangle t(m.vertices()[index1], m.vertices()[index2], m.vertices()[index3] );
        m.push_back(t);
    }

    return m;
}
#endif  //! OPENMEEG_CGAL_MESH_INCLUDE_H
