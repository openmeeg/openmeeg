// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <string>

#include "mesh.h"
#include "geometry.h"
#include "commandline.h"

using namespace OpenMEEG;

int
main(int argc,char* argv[]) {

    print_version(argv[0]);

    const CommandLine cmd(argc,argv,"Print Geometry information");
    const std::string& geom_filename = cmd.option("-g",std::string(),"Input .geom file");
    const std::string& cond_filename = cmd.option("-c",std::string(),"Input .cond file");

    if (cmd.help_mode())
        return 0;

    if (geom_filename=="") {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    Geometry geo;
    try {
        if (cond_filename!="") {
            geo.load(geom_filename,std::string(cond_filename));
        } else {
            geo.load(geom_filename);
        }
    } catch (maths::Exception& e) {
        std::cerr << e.what() << std::endl;
        return e.code();
    } catch (OpenMEEG::Exception& e) {
        std::cerr << e.what() << std::endl;
        return e.code();
    }

    if (!geo.selfCheck())
        return 1;

    std::cout << ".geom : OK" << std::endl;
    return 0;
}
