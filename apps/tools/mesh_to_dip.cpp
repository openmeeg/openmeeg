// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include "mesh.h"
#include "commandline.h"
#include "matrix.h"

using namespace OpenMEEG;

int
main(int argc,char* argv[]) {

    print_version(argv[0]);

    const CommandLine cmd(argc,argv,"Convert mesh file to a dipole file");
    const std::string& input_filename  = cmd.option("-i",std::string(),"Input Mesh");
    const std::string& output_filename = cmd.option("-o",std::string(),"Output .dip file");

    if (cmd.help_mode())
        return 0;

    if (input_filename=="" || output_filename=="") {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    Mesh m(input_filename,false);

    Matrix mat(m.vertices().size(),6);

    unsigned i = 0;
    for (const auto& vertex : m.vertices()) {
        mat(i,0) = vertex->x();
        mat(i,1) = vertex->y();
        mat(i,2) = vertex->z();
        const Normal& n = m.normal(*vertex);
        mat(i,3)   = n.x();
        mat(i,4)   = n.y();
        mat(i++,5) = n.z();
    }

    mat.save(output_filename);

    return 0;
}
