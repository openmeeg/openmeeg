// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the GPL open source license.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <OM_CGAL.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>

namespace OpenMEEG::OMCGAL {

    namespace Details {

        class LevelSet {
        public:

            LevelSet(const std::string& filename,const double levelset,const bool inside,const double outside_value):
                level(levelset),positive_inside(inside),outside(outside_value)
            {
                image.read(filename.c_str());
            }

            double operator()(const Point3& p) const {
                const double v = image.value(p.x(),p.y(),p.z());

                const double d = v-level;

                return (positive_inside) ? d : -d;
            }

            int xdim() const { return image.xdim(); }
            int ydim() const { return image.ydim(); }
            int zdim() const { return image.zdim(); }

        private:

            using Image = CGAL::Image_3;

            Image  image;
            double level;
            bool   positive_inside;
            double outside;
        };
    }

    Mesh mesh_from_levelset_image(const std::string& input_filename,const double levelset_value,const bool positive_inside,const double radius_bound,const double distance_bound) {

        namespace params = CGAL::parameters;

        using MeshDomain    = CGAL::Labeled_mesh_domain_3<Kernel>;
        using Triangulation = CGAL::Mesh_triangulation_3<MeshDomain>::type;
        using Criteria      = CGAL::Mesh_criteria_3<Triangulation>;
        using C3t3          = CGAL::Mesh_complex_3_in_triangulation_3<Triangulation>;

        const Details::LevelSet image(input_filename,levelset_value,positive_inside,1.0);

        std::cout << "Input INR image:\n dimension: " << image.xdim() << "x"<< image.ydim() << "x"<< image.zdim() << "\n Positive values are " << (positive_inside?"Inside":"Outside") << std::endl;

        // Mesh criteria

        // Carefully choose bounding sphere: the center must be inside the
        // surface defined by 'image' and the radius must be high enough so that
        // the sphere actually bounds the whole image.

        const Point3  bounding_sphere_center(image.xdim()/2.,image.ydim()/2.,image.zdim()/2.);
        const FT      bounding_sphere_squared_radius = image.xdim()*image.ydim()*2.;
        const Sphere3 bounding_sphere(bounding_sphere_center,bounding_sphere_squared_radius);

        // Domain: definition of the surface, with 1e-3 as relative precision.

        const MeshDomain& domain = MeshDomain::create_implicit_mesh_domain(image,bounding_sphere,params::relative_error_bound(1e-3));

        Criteria crit(params::facet_angle(30).facet_size(radius_bound).facet_distance(distance_bound));

        // Meshing

        const C3t3& c3t3 = CGAL::make_mesh_3<C3t3>(domain,crit,params::no_exude().no_perturb());

        return CGAL_to_OM(c3t3);
    }
}
