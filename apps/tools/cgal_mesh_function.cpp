// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <mesh.h>
#include "options.h"
#include <cgal_lib.h>

using namespace OpenMEEG;

int main(int argc, char **argv) {
    command_usage("Create a BEM mesh from either an implicit function: sphere, hemisphere, ...:");
    const double sphere_radius     = command_option("-r", 0.0, "radius of the sphere");
    const double hemisphere_radius = command_option("-hr",0.0, "radius of the hemisphere");
    const double radius_bound      = command_option("-fs",1e-1,"facet radius bound of elements");
    const double distance_bound    = command_option("-fd",1e-1,"facet distance bound to the input surface");
    // const unsigned init_points  = command_option("-ip", 10, "initial number of points (for the hemisphere)");
    const char * output_filename   = command_option("-o",nullptr,"Output Mesh");

    if (command_option("-h",nullptr,nullptr))
        return 0;

    if (output_filename==nullptr) {
        std::cerr << "Set an output filename" << std::endl;
        return 0;
    }

    Mesh m_out = cgal_mesh_function(sphere_radius,hemisphere_radius,radius_bound,distance_bound);
    m_out.save(output_filename);
    m_out.info();

    return 0;
}
