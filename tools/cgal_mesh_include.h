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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_3/Robust_intersection_traits_3.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/make_surface_mesh.h>

using namespace OpenMEEG;

// Domain 
// (we use exact intersection computation with Robust_intersection_traits_3)
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Mesh_3::Robust_intersection_traits_3<K> Geom_traits;
typedef Geom_traits::Point_3 Point_3;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Gray_level_image_3<Geom_traits::FT, Point_3> Gray_level_image;

typedef CGAL::Simple_cartesian<double> Simple_cartesian_kernel;
typedef CGAL::Implicit_surface_3<Geom_traits, Gray_level_image> Surface_3;

// Triangulation
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;
typedef Tr::Vertex_handle Vertex_handle;

// from polygons
namespace CGAL {
   template <class Polyhedron, class Kernel, class Dummy_kernel> class AABB_polyhedral_oracle;
}

typedef CGAL::AABB_polyhedral_oracle<Polyhedron,K,Simple_cartesian_kernel> Input_surface;

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

namespace CGAL {
    // class dealing a surface as an input for the surface mesher (or building an oracle from the input surface)
    template <class Polyhedron, class Kernel, class Dummy_kernel>
    class AABB_polyhedral_oracle  {
    public:
        typedef typename Kernel::FT FT;
        typedef typename Kernel::Ray_3 Ray_3;
        typedef typename Kernel::Line_3 Line_3;
        typedef typename Kernel::Point_3 Point_3;
        typedef typename Kernel::Segment_3 Segment_3;

        typedef AABB_polyhedral_oracle<Polyhedron,Kernel,Dummy_kernel> Self;
        typedef Self Surface_mesher_traits_3;
        typedef Point_3 Intersection_point;
        typedef Self Surface_3;

        // AABB tree
        typedef AABB_const_polyhedron_triangle_primitive<Kernel, Polyhedron> AABB_primitive;
        typedef class AABB_traits<Kernel,AABB_primitive> AABB_traits;
        typedef AABB_tree<AABB_traits> Tree;
        typedef typename AABB_traits::Bounding_box Bounding_box;

        typedef boost::shared_ptr<Tree> Tree_shared_ptr;
        Tree_shared_ptr m_pTree;

        Tree* tree() const { return m_pTree.get(); }

        // Surface constructor
        AABB_polyhedral_oracle(const Polyhedron& poly) : m_pTree(Tree_shared_ptr(new Tree(poly.facets_begin(), poly.facets_end()))) { }

        AABB_polyhedral_oracle(const AABB_polyhedral_oracle& oracle) : m_pTree(oracle.m_pTree) { }

        class Intersect_3;
        friend class Intersect_3;
        class Intersect_3 {
            typedef boost::optional<typename Tree::Object_and_primitive_id> 
                AABB_intersection;

            const Self& self;
            public:
            Intersect_3(const Self& self_) : self(self_) { }

            Object operator()(const Surface_3& surface, const Segment_3& segment) const
            {
                AABB_intersection intersection = surface.tree()->any_intersection(segment);

                if ( intersection )
                    return intersection->first;
                else
                    return Object();
            }

            Object operator()(const Surface_3& surface, const Ray_3& ray) const
            {
                AABB_intersection intersection = surface.tree()->any_intersection(ray);

                if ( intersection )
                    return intersection->first;
                else
                    return Object();
            }

            Object operator()(const Surface_3& surface, const Line_3& line) const
            {
                AABB_intersection intersection = surface.tree()->any_intersection(line);

                if ( intersection )
                    return intersection->first;
                else
                    return Object();
            }
        };
        Intersect_3 intersect_3_object() const
        {
            return Intersect_3(*this);
        }

        class Construct_initial_points;
        friend class Construct_initial_points;
        class Construct_initial_points
        {
            const Self& self;
            public:
            Construct_initial_points(const Self& self_) : self(self_)
            {
            }

            template <typename OutputIteratorPoints>
                OutputIteratorPoints operator() (const Surface_3, OutputIteratorPoints out, int ) const
                {
                    return out;
                }
        };
        Construct_initial_points construct_initial_points_object() const
        {
            return Construct_initial_points(*this);
        }

        template <class P>
            bool is_in_volume(const Surface_3& surface, const P& p)
            {
                const Bounding_box bbox = surface.tree()->root_bbox();
                if(p.x() < bbox.xmin() || p.x() > bbox.xmax())
                    return false;
                if(p.y() < bbox.ymin() || p.y() > bbox.ymax())
                    return false;
                if(p.z() < bbox.zmin() || p.z() > bbox.zmax())
                    return false;

                const double diameter = max_bbox_length(bbox) * 2;

                typename CGAL::Random_points_on_sphere_3<Point_3> random_point(FT(1));
                typename Kernel::Construct_vector_3 vector =
                    Kernel().construct_vector_3_object();
                typename Kernel::Construct_segment_3 segment =
                    Kernel().construct_segment_3_object();
                typename Kernel::Construct_translated_point_3 translate =
                    Kernel().construct_translated_point_3_object();
                typename Kernel::Construct_scaled_vector_3 scale =
                    Kernel().construct_scaled_vector_3_object();

                const Segment_3 seg =
                    segment(p, translate(p,
                                scale(vector(ORIGIN,
                                        *random_point),
                                    diameter)));
                return (surface.tree()->number_of_intersections(seg) % 2) == 1;
            }
    private:
        double max_bbox_length(const Bounding_box& bbox) const
        {
            return (std::max)(bbox.xmax()-bbox.xmin(),
                    (std::max)(bbox.ymax()-bbox.ymin(),
                        bbox.zmax()-bbox.zmin()));
        }
    }; // end class AABB_polyhedral_oracle
} // end namespace CGAL

Mesh CGAL_to_OM(C2t3 c2t3) {
    Mesh m(c2t3.triangulation().number_of_vertices(), c2t3.number_of_facets());
    bool success;

    std::map<Vertex_handle, unsigned> V;
    unsigned inum=0;
    for(Tr::Finite_vertices_iterator vit = c2t3.triangulation().finite_vertices_begin(), end = c2t3.triangulation().finite_vertices_end(); vit != end; ++vit)
    {
        const Point_3& p = vit->point();
        Vect3 pts(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
        m[inum] = pts;
        V[vit] = inum++;
    }

    Tr::Finite_facets_iterator fit = c2t3.triangulation().finite_facets_begin();
    C2t3::Facets oriented_set;
    std::stack<C2t3::Facet> stack;

    Tr::size_type number_of_facets = c2t3.number_of_facets();

    while (oriented_set.size() != number_of_facets) 
    {
        while ( fit->first->is_facet_on_surface(fit->second) == false || oriented_set.find(*fit) != oriented_set.end() || oriented_set.find(c2t3.opposite_facet(*fit)) != oriented_set.end() ) {
            ++fit;
        }
        oriented_set.insert(*fit);
        stack.push(*fit);
        while(! stack.empty() ) {
            C2t3::Facet f = stack.top();
            stack.pop();
            for(int ih = 0 ; ih < 3 ; ++ih) {
                const int i1  = c2t3.triangulation().vertex_triple_index(f.second, c2t3.triangulation().cw(ih));
                const int i2  = c2t3.triangulation().vertex_triple_index(f.second, c2t3.triangulation().ccw(ih));
                const C2t3::Face_status face_status = c2t3.face_status(C2t3::Edge(f.first, i1, i2));
                if (face_status == C2t3::REGULAR) {
                    C2t3::Facet fn = c2t3.neighbor(f, ih);
                    if (oriented_set.find(fn) == oriented_set.end()) {
                        if(oriented_set.find(c2t3.opposite_facet(fn)) == oriented_set.end())
                        {
                            oriented_set.insert(fn);
                            stack.push(fn);
                        }
                        else {
                            success = false; // non-orientable
                        }
                    }
                }
                else if(face_status != C2t3::BOUNDARY) {
                    success = false; // non manifold, thus non-orientable
                }
            } // end "for each neighbor of f"
        } // end "stack non empty"
    } // end "oriented_set not full"
    inum = 0;
    
    for(C2t3::Facets::const_iterator fcit = oriented_set.begin(); fcit != oriented_set.end(); ++fcit) {
      const Tr::Cell_handle cell = fcit->first;
      const int& index = fcit->second;
      const int index1 = V[cell->vertex(c2t3.triangulation().vertex_triple_index(index, 0))];
      const int index2 = V[cell->vertex(c2t3.triangulation().vertex_triple_index(index, 1))];
      const int index3 = V[cell->vertex(c2t3.triangulation().vertex_triple_index(index, 2))];
      Triangle t(index1, index2, index3, Vect3(0.,0.,0.));
      m.triangle(inum++) = t;
    }

    return m;
}
#endif  //! OPENMEEG_CGAL_MESH_INCLUDE_H
