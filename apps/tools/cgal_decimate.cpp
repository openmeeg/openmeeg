// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

// Mesh simplification adapted from Kai Dang 2015 @ Inria

#include <mesh.h>
#include <cgal_lib.h>

#include "commandline.h"

using namespace OpenMEEG;

int main(int argc,char **argv) {

    const CommandLine cmd(argc,argv,"Decimate a mesh:");
    const std::string input_filename  = cmd.option("-i",std::string(),"Input image or mesh");
    const unsigned    nb_points       = cmd.option("-n",1000,"desired output number of vertices");
    const std::string output_filename = cmd.option("-o",std::string(),"Output Mesh");

    if (cmd.help_mode())
        return 0;

    const Mesh mesh(input_filename,false);
    std::cout << "Input surface:\n nb of points: " << mesh.vertices().size() << "\t nb of triangles:\t" << mesh.triangles().size() << std::endl;

    const Mesh& decimated_mesh = cgal_decimate(mesh,nb_points);
    decimated_mesh.info();
    decimated_mesh.save(output_filename);

    return 0;
}
