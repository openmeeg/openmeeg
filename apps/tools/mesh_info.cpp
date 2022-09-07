// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include "mesh.h"
#include "logger.h"
#include "commandline.h"

using namespace OpenMEEG;

int
main(int argc,char* argv[]) {

    print_version(argv[0]);

    const CommandLine cmd(argc,argv,"Get info about a Mesh");
    const std::string& input_filename = cmd.option("-i",std::string(),"Input Mesh");

    if (cmd.help_mode())
        return 0;

    if (input_filename=="") {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    Mesh m(input_filename);

    if (m.has_self_intersection())
        log_stream(WARNING) << "Mesh is self intersecting !" << std::endl;

    // for closed mesh

    if (!m.has_correct_orientation()) {
        log_stream(WARNING) << "Mesh is not well-oriented (valid for closed mesh) !" << std::endl;
        return 1;
    }

    //  For closed meshes E = 3*F/2
    //  For a simple closed surface, V-E+F=2.
    //  This the test for a closed mesh is V-F/2=2 or 2*V-F=4.

    if (2*m.vertices().size()-m.triangles().size()==4) {
        std::cout << "Mesh orientation correct (valid for closed mesh)." << std::endl;
    } else {
        std::cout << "Mesh local orientation correct." << std::endl;
    }

    return 0;
}
