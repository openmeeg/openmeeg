// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <filesystem>

#include "geometry.h"
#include "filenames.h"
#include "matrix.h"
#include "commandline.h"

using namespace OpenMEEG;

int
main(int argc,char* argv[]) {

    print_version(argv[0]);

    const CommandLine cmd(argc,argv,"Convert a single VTK/VTP into meshes.");
    const std::string& input  = cmd.option("-i",std::string(),"Input VTK/VTP file");
    const std::string& output = cmd.option("-o",std::string(),"Output mesh base name");

    if (cmd.help_mode())
        return 0;

    if (input=="" || output=="") {
        std::cout << "Missing arguments, try the -h option" << std::endl;
        return 1;
    }

    Geometry geom(input);

    const std::string  basename  = fs::path(output).stem();
    const std::string& extension = getFilenameExtension(output);

    for (const auto& mesh : geom.meshes())
        mesh.save(basename+mesh.name()+'.'+extension);

    geom.save(basename+".geom");

    return 0;
}
