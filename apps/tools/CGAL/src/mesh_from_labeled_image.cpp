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

    Mesh mesh_from_labeled_image(const std::string& input_filename,const double radius_bound,const double distance_bound) {

        namespace params = CGAL::parameters;

        using Image         = CGAL::Image_3;
        using MeshDomain    = CGAL::Labeled_mesh_domain_3<Kernel>;
        using Triangulation = CGAL::Mesh_triangulation_3<MeshDomain>::type;
        using Criteria      = CGAL::Mesh_criteria_3<Triangulation>;
        using C3t3          = CGAL::Mesh_complex_3_in_triangulation_3<Triangulation>;

        // Mesh criteria

        Image image;
        image.read(input_filename.c_str());

        std::cout << "Input image:\n dimension: " << image.xdim() << "x" << image.ydim() << "x" << image.zdim() << std::endl;

        // Domain and criteria.

        MeshDomain domain = MeshDomain::create_labeled_image_mesh_domain(image);
        Criteria crit(params::facet_angle(30).facet_size(radius_bound).facet_distance(distance_bound));

        // Meshing.

        const C3t3& c3t3 = CGAL::make_mesh_3<C3t3>(domain,crit,params::no_exude().no_perturb());

        return CGAL_to_OM(c3t3);
    }
}
