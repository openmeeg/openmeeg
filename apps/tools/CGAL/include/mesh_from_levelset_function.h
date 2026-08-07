// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the GPL open source license.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

// for verbosity
#define CGAL_MESH_3_VERBOSE

#include <string>

#include <mesh.h>

#include "OM_CGAL.h"
#include "Conversions.h"

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>

namespace OpenMEEG::OMCGAL {

    //  Implicit Sphere and HemiSphere functions.

    class SphereFunction {
    public:

        SphereFunction(const FT sqradius): squared_radius(sqradius) { }

        FT operator()(const Point3& p) const {
            using CGAL::square;
            return square(p.x())+square(p.y())+square(p.z())-squared_radius;
        }

    private:

        FT squared_radius;
    };

    class HemiSphereFunction {
    public:

        HemiSphereFunction(const FT sqradius): squared_radius(sqradius) { }

        FT operator()(const Point3& p) const {
            using CGAL::square;
            const double x2 = square(p.x());
            const double y2 = square(p.y());
            const double z2 = square(p.z());
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

    /// Create a mesh from an implicit function.

    template <typename Function,typename BoundingSphere>
    Mesh mesh_from_levelset_function(const Function& f,const BoundingSphere& bound,const double radius_bound,const double distance_bound) {

        namespace params = CGAL::parameters;

        using MeshDomain    = CGAL::Labeled_mesh_domain_3<Kernel>;
        using Triangulation = CGAL::Mesh_triangulation_3<MeshDomain>::type;
        using C3t3S         = CGAL::Mesh_complex_3_in_triangulation_3<Triangulation>;
        using Criteria      = CGAL::Mesh_criteria_3<Triangulation>;

        // Defining the domain,

        const MeshDomain& domain = MeshDomain::create_implicit_mesh_domain(f,bound,params::relative_error_bound(1e-6));

        // Mesh criteria

        Criteria crit(params::facet_angle(30).facet_size(radius_bound).facet_distance(distance_bound));

        // meshing domain

        const C3t3S& c3t3 = CGAL::make_mesh_3<C3t3S>(domain,crit,params::no_exude().no_perturb());

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
}
