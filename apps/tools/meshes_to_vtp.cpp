// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <geometry.h>
#include <matrix.h>
#include <commandline.h>

using namespace OpenMEEG;

inline std::string
param(const char* base,const unsigned n) {
    std::stringstream res;
    res << base << n;
    return res.str();
}

int main( int argc, char **argv) {

    print_version(argv[0]);

    const CommandLine cmd(argc,argv,"Convert meshes OR geometry into a single VTK/VTP file.");
    const std::string& geom_filename = cmd.option("-g",std::string(),"Input geometry");

    Geometry::MeshList meshes;
    for (unsigned i=0; i<7; ++i) {
        const std::string& mesh_path_option_name = param("-i",i+1);
        const std::string& mesh_name_option_name = param("-n",i+1);
        const std::string& mesh_name_default_val = param("",i+1);
        const std::string& path = cmd.option(mesh_path_option_name,std::string(),        "Input mesh");
        const std::string& name = cmd.option(mesh_name_option_name,mesh_name_default_val,"Mesh name");
        meshes.push_back({ name, path });
    }

    const std::string& output = cmd.option("-o",std::string(),"Output VTP file");

    if (cmd.help_mode())
        return 0;

    if ((meshes.size()==0 && geom_filename=="") || output=="") {
        std::cout << "Missing arguments, try the -h option" << std::endl;
        return 1;
    }
    
    Geometry geom;

    if (geom_filename!="") {
        geom.load(geom_filename);
    } else {
        geom.import(meshes);
    }

    geom.info();
    geom.save(output);

    return 0;
}
