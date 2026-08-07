// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the GPL open source license.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <mesh.h>
#include <OM_CGAL.h>

#include "commandline.h"

using namespace OpenMEEG;

int main(int argc,char **argv) {

    using CGAL::ORIGIN;
    using CGAL::square;

    using FT      = OpenMEEG::OMCGAL::FT;
    using Point3  = OpenMEEG::OMCGAL::Point3;
    using Sphere3 = OpenMEEG::OMCGAL::Sphere3;

    CommandLine cmd(argc,argv,"Create a BEM mesh from either an implicit function: sphere, hemisphere, ...:");
    const double       sphere_radius     = cmd.option("-r", 0.0, "radius of the sphere");
    const double       hemisphere_radius = cmd.option("-hr",0.0, "radius of the hemisphere");
    const double       radius_bound      = cmd.option("-fs",1e-1,"facet radius bound of elements");
    const double       distance_bound    = cmd.option("-fd",1e-1,"facet distance bound to the input surface");
    // const unsigned init_points  = cmd.option("-ip", 10, "initial number of points (for the hemisphere)");
    const std::string& filename          = cmd.option("-o",std::string(),"Output Mesh");

    if (cmd.help_mode())
        return 0;

    if (filename.empty()) {
        std::cerr << argv[0] << " needs an output filename" << std::endl;
        return 1;
    }

    if (sphere_radius!=0.0 && hemisphere_radius!=0.0) {
        std::cerr << "Only one of -r and -hr can be given" << std::endl;
        return 2;
    }

    if (sphere_radius==0.0 && hemisphere_radius==0.0) {
        std::cerr << "At least one of -r and -hr must be given" << std::endl;
        return 3;
    }

    auto make_mesh = [&radius_bound,&distance_bound,&filename]<typename Function>(const Function& function,const Sphere3& bound) {
        const Mesh& mesh = OpenMEEG::OMCGAL::mesh_from_levelset_function(function,bound,radius_bound,distance_bound);
        mesh.save(filename);
        mesh.info();
    };

    if (sphere_radius!=0) { // Sphere case.

        using SphereFunction = OpenMEEG::OMCGAL::SphereFunction;

        const SphereFunction spherefunction(square(sphere_radius));
        const Sphere3        bounding_sphere(ORIGIN,square(FT(1.1)*sphere_radius));

        make_mesh(spherefunction,bounding_sphere);

    } else { // Hemisphere case.

        using HemiSphereFunction = OpenMEEG::OMCGAL::HemiSphereFunction;

        const HemiSphereFunction hemispherefunction(hemisphere_radius);
        const Sphere3            bounding_sphere(Point3(0.0,0.0,hemisphere_radius/2.0),square(FT(1.1)*sphere_radius));

        make_mesh(hemispherefunction,bounding_sphere);
    }

    return 0;
}
