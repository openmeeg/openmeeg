// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the GPL open source license.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <string>

#include <mesh.h>
#include <OM_CGAL.h>

#include "commandline.h"

using namespace OpenMEEG;

int main(int argc,char **argv) {

    CommandLine cmd(argc,argv,"Re-mesh a mesh:");
    const std::string& input_filename  = cmd.option("-i",std::string(),"Input image or mesh");
    const double       radius_bound    = cmd.option("-fs",1e-1,"facet radius bound of elements");
    const double       distance_bound  = cmd.option("-fd",1e-1,"facet distance bound to the input surface");
    const std::string& output_filename = cmd.option("-o",std::string(),"Output Mesh");
    const std::string& sizing_field    = cmd.option("-field",std::string(),"(OPTIONAL) definition of the space to be refined 3 times finer (a matrix file: with either: \"x y z nx ny nz\" per line to define planes (by intersection of domains), or \"x y z r\" to define spheres (by union of domains).)");

    if (cmd.help_mode())
        return 0;

    if (output_filename.empty()) {
        std::cerr << argv[0] << " needs an output filename." << std::endl;
        return 1;
    }

    // Mesh input

    const Mesh mesh(input_filename,false);
    std::cout << "Input surface:\n nb of points: " << mesh.vertices().size() << "\t nb of triangles:\t" << mesh.triangles().size() << std::endl;

    const Mesh& refined_mesh = OpenMEEG::OMCGAL::refine_mesh(mesh,radius_bound,distance_bound,sizing_field);
    refined_mesh.save(output_filename);
    refined_mesh.info();

    return 0;
}
