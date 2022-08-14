// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include "mesh.h"
#include "commandline.h"

using namespace OpenMEEG;

int
main( int argc,char* argv[]) {

    print_version(argv[0]);

    const CommandLine cmd(argc,argv,"Concat 2 mesh and save the result");
    
    const std::string& input_filename1 = cmd.option("-i1",std::string(),"Input Mesh 1");
    const std::string& input_filename2 = cmd.option("-i2",std::string(),"Input Mesh 2");
    const std::string& output_filename = cmd.option("-o", std::string(),"Output Mesh");

    if (cmd.help_mode())
        return 0;

    if (argc<2) {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    Mesh m1(input_filename1);

    Mesh m2(input_filename2);
    Mesh m3;

    m3.merge(m1, m2);
    m3.save(output_filename);

    return 0;
}
