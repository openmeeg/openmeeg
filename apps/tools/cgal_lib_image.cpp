// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <cgal_lib.h>

#include <CGAL/Labeled_image_mesh_domain_3.h>

namespace OpenMEEG {

    Mesh cgal_mesh_3Dlabeled_image(const char* input_filename,const double radius_bound,const double distance_bound) {

        using Image     = CGAL::Image_3;
        using Domain    = CGAL::Labeled_image_mesh_domain_3<Image>;
        using Criterion = CGAL::Mesh_criteria_3<CGAL::Mesh_triangulation_3<Domain>::type>;

        // Mesh criteria

        Image image;
        image.read(input_filename);

        std::cout << "Input image:\n dimension: " << image.xdim() << "x" << image.ydim() << "x" << image.zdim() << std::endl;

        // Domain and criteria.

        Domain    domain(image,0);
        Criterion crit(facet_angle=30,facet_size=radius_bound,facet_distance=distance_bound);

        // Meshing.

        return CGAL_to_OM(CGAL::make_mesh_3<C3t3>(domain,crit,CGAL::parameters::no_exude(),CGAL::parameters::no_perturb()));
    }

    Mesh cgal_mesh_3Dlevelset_image(const char* input_filename,const double levelset_value,const bool positive_inside,const double radius_bound,const double distance_bound) {

        using Image     = CGAL::Gray_level_image_3<K::FT,Point_3>; // Grey level image.
        using Domain    = CGAL::Implicit_mesh_domain_3<Image,K>;
        using Criterion = CGAL::Mesh_criteria_3<CGAL::Mesh_triangulation_3<Domain>::type>;

        // Mesh criteria

        double value_outside = 1.;
        Gray_level_image image(input_filename,levelset_value,positive_inside,value_outside);

        std::cout << "Input INR image:\n dimension: " << image.xdim() << "x"<< image.ydim() << "x"<< image.zdim() << "\n Positive values are " << (positive_inside?"Inside":"Outside") << std::endl;

        // Carefully choose bounding sphere: the center must be inside the
        // surface defined by 'image' and the radius must be high enough so that
        // the sphere actually bounds the whole image.

        Point_3 bounding_sphere_center(image.xdim()/2.,image.ydim()/2.,image.zdim()/2.);
        K::FT bounding_sphere_squared_radius = image.xdim()*image.ydim()*2.;
        K::Sphere_3 bounding_sphere(bounding_sphere_center,bounding_sphere_squared_radius);

        // Domain: definition of the surface, with 1e-8 as relative precision.

        Domain    domain(image,bounding_sphere,1e-3);
        Criterion crit(facet_angle=30,facet_size=radius_bound,facet_distance=distance_bound);

        // Meshing

        return CGAL_to_OM(CGAL::make_mesh_3<C3t3L>(domain,crit,CGAL::parameters::no_exude(),CGAL::parameters::no_perturb()));
    }
}
