// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the GPL open source license.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <string>

#include <OM_CGAL.h>

#include "commandline.h"

using namespace OpenMEEG;

int main(int argc,char **argv) {

    CommandLine cmd(argc,argv,"Create a BEM mesh from a 3D levelset image:");
    const std::string& input_filename  = cmd.option("-i",std::string(),"Input image");
    const double       levelset_value  = cmd.option("-v",0.0,"Levelset value");
    const bool         inout           = cmd.option("-inout",false,"Inside out the image ?");
    const double       radius_bound    = cmd.option("-fs",1e-1,"facet radius bound of elements");
    const double       distance_bound  = cmd.option("-fd",1e-1,"facet distance bound to the input surface");
    const std::string& output_filename = cmd.option("-o",std::string(),"Output Mesh");

    if (cmd.help_mode())
        return 0;

    if (output_filename.empty()) {
        std::cerr << argv[0] << " needs an output filename." << std::endl;
        return 1;
    }

    // Output

    const Mesh& mesh = OpenMEEG::OMCGAL::mesh_from_levelset_image(input_filename,levelset_value,inout,radius_bound,distance_bound);

    mesh.save(output_filename);
    mesh.info();

    return 0;
}
