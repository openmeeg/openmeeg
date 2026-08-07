// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the GPL open source license.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include "OM_CGAL.h"
#include "Conversions.h"

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>

// Simplification function

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>

namespace OpenMEEG::OMCGAL {

    using PolyhedralDomain = CGAL::Polyhedral_mesh_domain_3<SurfaceMesh,Kernel>;

    namespace Details {

        // Sizing field

        struct Planes {

            Planes(const OpenMEEG::Matrix& M,const double field_size): mat(M),fs(field_size) { }

            FT operator()(const Point3& p,const int,const PolyhedralDomain::Index&) const {
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

            FT operator()(const Point3& p,const int,const PolyhedralDomain::Index&) const {
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
    }

    /// Refine a mesh.

    Mesh refine_mesh(const Mesh& mesh,const double radius_bound,const double distance_bound,const std::string& sizing_field) {

        namespace params = CGAL::parameters;

        using Triangulation = CGAL::Mesh_triangulation_3<PolyhedralDomain>::type;
        using C3t3          = CGAL::Mesh_complex_3_in_triangulation_3<Triangulation>;
        using Criteria      = CGAL::Mesh_criteria_3<Triangulation>;

        // Create input mesh and domain.

        SurfaceMesh      cgal_mesh = OM_to_CGAL_SurfaceMesh(mesh);
        PolyhedralDomain domain(cgal_mesh);

        // Create the refined mesh.

        C3t3 c3t3;

        if (sizing_field.empty()) {

            const Criteria criteria(params::facet_angle(30).facet_size(radius_bound).facet_distance(distance_bound));
            c3t3 = CGAL::make_mesh_3<C3t3>(domain,criteria,params::no_exude().no_perturb());

        } else {

            Matrix field(sizing_field);
            if (field.ncol()==6) {
                const Details::Planes planes(field,radius_bound);
                const Criteria criteria(params::facet_angle(30).facet_size(planes).facet_distance(distance_bound));
                c3t3 = CGAL::make_mesh_3<C3t3>(domain,criteria,params::no_exude().no_perturb());
            } else if (field.ncol()==4) {
                const Details::Spheres spheres(field,radius_bound);
                const Criteria criteria(params::facet_angle(30).facet_size(spheres).facet_distance(distance_bound));
                c3t3 = CGAL::make_mesh_3<C3t3>(domain,criteria,params::no_exude().no_perturb());
            } else {
                std::cerr << "Error: file should contain either 4 or 6 columns" << std::endl;
            }
        }

        return CGAL_to_OM(c3t3);
    }
}
