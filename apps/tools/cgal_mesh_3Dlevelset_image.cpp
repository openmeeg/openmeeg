// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <options.h>
#include <cgal_lib.h>

using namespace OpenMEEG;

int main(int argc, char **argv) {
    command_usage("Create a BEM mesh from a 3D levelset image:");
    const char * input_filename  = command_option("-i",(const char *) NULL,"Input image");
    const double levelset_value  = command_option("-v", 0.,"Levelset value");
    const bool   inout           = command_option("-inout",false,"Inside out the image ?");
    const double radius_bound    = command_option("-fs",1e-1,"facet radius bound of elements");
    const double distance_bound  = command_option("-fd",1e-1,"facet distance bound to the input surface");
    const char * output_filename = command_option("-o",(const char *) NULL,"Output Mesh");

    if ( command_option("-h",(const char *)0,0) ) {
        return 0;
    }
    if ( output_filename == NULL ) {
        std::cerr << "Set an output filename" << std::endl;
        return 0;
    }

    // Output
    Mesh m_out = cgal_mesh_3Dlevelset_image(input_filename, levelset_value, inout, radius_bound, distance_bound);
    m_out.save(output_filename);
    m_out.info();
    return 0;
}
