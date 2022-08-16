// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <iostream>
#include "matrix.h"
#include "symmatrix.h"
#include "vector.h"
#include "commandline.h"

using namespace OpenMEEG;

int
main(int argc,char** argv) {

    print_version(argv[0]);

    const CommandLine cmd(argc,argv,"Convert squids in text file to a vtk file for vizualisation");
    const std::string& input_filename  = cmd.option("-i",std::string(),"Squids positions in original coordinate system");
    const std::string& output_filename = cmd.option("-o",std::string(),"Squids positions with orientations in vtk format");
    if (cmd.help_mode())
        return 0;

    if(input_filename=="" || output_filename=="") {
        std::cerr << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    Matrix squids(input_filename);

    std::ofstream ofs(output_filename);
    if (!ofs) {
        std::cerr << "Cannot open file " << output_filename << " for writing." << std::endl;
        return 1;
    }

    ofs <<"# vtk DataFile Version 3.0" << std::endl
        << "vtk output" << std::endl
        << "ASCII" << std::endl
        << "DATASET POLYDATA" << std::endl
        << "POINTS " << squids.nlin() << " float" << std::endl;

    for (unsigned int i=0; i<squids.nlin(); ++i)
        ofs << squids(i,0) << ' ' << squids(i,1) << ' ' << squids(i,2) << std::endl;

    ofs << "POINT_DATA " << squids.nlin() << std::endl
        << "NORMALS normals float" << std::endl;

    for (unsigned int i=0; i<squids.nlin(); ++i)
        ofs << squids(i,3) << ' ' << squids(i,4) << ' ' << squids(i,5) << std::endl;

    return 0;
}
