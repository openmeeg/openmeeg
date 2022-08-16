// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <string>

#include <matrix.h>
#include <danielsson.h>
#include <vector.h>
#include <sensors.h>
#include <commandline.h>

using namespace OpenMEEG;

int
main(int argc,char* argv[]) {

    print_version(argv[0]);

    const CommandLine cmd(argc,argv,"Project the sensors onto the given mesh:");
    const std::string& sensors_filename = cmd.option("-i",std::string(),"Sensors positions");
    const std::string& mesh_filename    = cmd.option("-m",std::string(),"Mesh on which to project the sensors");
    const std::string& output_filename  = cmd.option("-o",std::string(),"Output sensors positions");

    if (cmd.help_mode())
        return 0;

    if (sensors_filename=="" || mesh_filename=="" || output_filename=="") {
        std::cout << "Missing arguments, try the -h option" << std::endl;
        return 1;
    }

    // Read the file containing the positions of the EEG patches

    Sensors sensors(sensors_filename);

    Mesh mesh(mesh_filename);
    Interface interface;
    interface.oriented_meshes().push_back(OrientedMesh(mesh,OrientedMesh::Normal)); // one mesh per interface, (well oriented)

    Matrix output(sensors.getNumberOfPositions(), 3);

    const size_t nb_positions = sensors.getNumberOfPositions();
    for (size_t i=0; i<nb_positions; ++i) {
        const Vector position = sensors.getPosition(i);
        Vect3 current_position;
        for (unsigned k=0; k<3; ++k)
            current_position(k) = position(k);
        Vect3 alphas;
        const auto& res = dist_point_interface(current_position,interface,alphas);
        const Triangle& triangle = std::get<1>(res); // closest triangle
        current_position = alphas(0)*triangle.vertex(0)+alphas(1)*triangle.vertex(1)+alphas(2)*triangle.vertex(2);
        for (unsigned k=0; k<3; ++k)
            output(i,k) = current_position(k);
    }

    sensors.getPositions() = output;
    sensors.save(output_filename);

    return 0;
}
