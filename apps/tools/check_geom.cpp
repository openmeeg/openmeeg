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

    const CommandLine cmd(argc,argv,"Check mesh intersections in geometry file");
    const std::string& geom_filename = cmd.option("-g",std::string(),"Input .geom file");
    const std::string& mesh_filename = cmd.option("-m",std::string(),"Mesh file (ex: to test .geom with cortex mesh)");
    const std::string& dip_filename  = cmd.option("-d",std::string(),"The dipole .dip file (ex: to test .geom with cortical dipoles");
    const bool         verbose       = cmd.option("-v",false,        "Print verbose information about the geometry");

    if (cmd.help_mode())
        return 0;

    if (geom_filename=="") {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    Geometry g(geom_filename);

    if (!g.selfCheck())
        return 1;

    if (verbose) {
        std::cout << "Detailed information about the geom file :" << std::endl;
        g.info(true);
    }

    std::cout << ".geom : OK" << std::endl;
    if (mesh_filename!="") {
        Mesh m(mesh_filename);
        if (!g.check(m))
            return 1;
        std::cout << ".geom and mesh : OK" << std::endl;
    }

    if (dip_filename!="") {
        if (!g.is_nested()) {
            std::cerr << "Dipoles are only allowed when geometry is nested." << std::endl;
            return 1;
        }
        Matrix dipoles(dip_filename);
        if (!g.check_inner(dipoles))
            return 1;
        std::cout << ".geom and .dip dipoles : OK" << std::endl;
    }
    return 0;
}
