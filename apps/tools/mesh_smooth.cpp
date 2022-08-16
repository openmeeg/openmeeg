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
main(int argc,char* argv[]) {

    print_version(argv[0]);

    const CommandLine cmd(argc,argv,"Get info about a Mesh");
    const std::string& input_filename      = cmd.option("-i",std::string(),"Input Mesh");
    const std::string& output_filename     = cmd.option("-o",std::string(),"Output Mesh");
    const double       smoothing_intensity = cmd.option("-s",0.1,          "Smoothing Intensity");
    const size_t       niter               = cmd.option("-n",1000,         "Number of iterations");

    if (cmd.help_mode())
        return 0;

    if (input_filename=="" || output_filename=="") {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    Mesh m(input_filename);
    m.smooth(smoothing_intensity, niter);
    std::cout << "Smoothing done !" << std::endl;
    m.info();
    m.save(output_filename);

    return 0;
}
