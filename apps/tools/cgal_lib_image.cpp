// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <cgal_lib.h>

#include <CGAL/Labeled_image_mesh_domain_3.h>
#include <CGAL/Gray_level_image_3.h>

typedef CGAL::Labeled_image_mesh_domain_3<CGAL::Image_3,K> Mesh_domain;
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

typedef CGAL::Mesh_3::Robust_intersection_traits_3<K> Geom_traits;
typedef CGAL::Gray_level_image_3<Geom_traits::FT, Geom_traits::Point_3> Gray_level_image;
typedef CGAL::Implicit_mesh_domain_3<Gray_level_image, K> Mesh_domainL;
typedef CGAL::Mesh_triangulation_3<Mesh_domainL>::type TrL;
typedef CGAL::Mesh_complex_3_in_triangulation_3<TrL> C3t3L;
// Criteria
typedef CGAL::Mesh_criteria_3<TrL> Mesh_criteriaL;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

namespace OpenMEEG {
    Mesh cgal_mesh_3Dlabeled_image(const char* input_filename, double radius_bound, double distance_bound)
    {
        // Mesh criteria
        CGAL::Image_3 image;
        image.read(input_filename);
        std::cout << "Input image:\n dimension: " << image.xdim() << "x"<< image.ydim() << "x"<< image.zdim() << std::endl;

        // Domain
        Mesh_domain domain(image, 0);

        // Mesh criteria
        Mesh_criteria criteria(facet_angle=30, facet_size=radius_bound, facet_distance=distance_bound);

        // Meshing
        C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());

        // Output
        return CGAL_to_OM(c3t3);
    }

    Mesh cgal_mesh_3Dlevelset_image(const char* input_filename,  double levelset_value, bool positive_inside, double radius_bound, double distance_bound)
    {
        // Mesh criteria
        double value_outside  = 1.;
        Gray_level_image image(input_filename, levelset_value, positive_inside, value_outside);
        std::cout << "Input INR image:\n dimension: " << image.xdim() << "x"<< image.ydim() << "x"<< image.zdim() << "\n Positive values are " << (positive_inside?"Inside":"Outside") << std::endl;
        // Carefully choose bounding sphere: the center must be inside the
        // surface defined by 'image' and the radius must be high enough so that
        // the sphere actually bounds the whole image.
        Point_3 bounding_sphere_center(image.xdim()/2., image.ydim()/2., image.zdim()/2.);
        K::FT bounding_sphere_squared_radius = image.xdim()*image.ydim()*2.;
        K::Sphere_3 bounding_sphere(bounding_sphere_center, bounding_sphere_squared_radius);

        // Domain
        // definition of the surface, with 10^-8 as relative precision
        Mesh_domainL domain(image, bounding_sphere, 1e-3);

        // Mesh criteria
        Mesh_criteriaL criteria(facet_angle=30, facet_size=radius_bound, facet_distance=distance_bound);

        // Meshing
        C3t3L c3t3 = CGAL::make_mesh_3<C3t3L>(domain, criteria, no_exude(), no_perturb());

        // Output
        return CGAL_to_OM(c3t3);
    }
}
