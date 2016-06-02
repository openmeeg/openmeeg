/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

The OpenMEEG software is a C++ package for solving the forward/inverse
problems of electroencephalography and magnetoencephalography.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's authors,  the holders of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

#include <string>
#include <options.h>
#include <matrix.h>
#include <danielsson.h>
#include <vector.h>
#include <om_utils.h>
#include <sensors.h>

using namespace OpenMEEG;

int main( int argc, char** argv)
{
    print_version(argv[0]);

    command_usage("Project the sensors onto the given mesh:");
    const char *sensors_filename = command_option("-i",(const char *) NULL,"Sensors positions");
    const char *mesh_filename = command_option("-m",(const char *) NULL,"Mesh on which to project the sensors");
    const char *output_filename = command_option("-o",(const char *) NULL,"Output sensors positions");

    if (command_option("-h",(const char *)0,0)) return 0;

    if ( !sensors_filename || !mesh_filename || !output_filename) {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    // read the file containing the positions of the EEG patches
    Sensors sensors(sensors_filename);

    Mesh mesh(mesh_filename);
    Interface interface;
    interface.push_back(OrientedMesh(mesh, true)); // one mesh per interface, (well oriented)

    Matrix output(sensors.getNumberOfPositions(), 3);

    size_t nb_positions = sensors.getNumberOfPositions();

    for( size_t i = 0; i < nb_positions; ++i )
    {
        Vector position = sensors.getPosition(i);
        Vect3 current_position, alphas;
        for ( unsigned k = 0; k < 3; ++k) {
            current_position(k) = position(k);
        }
        Triangle triangle; // closest triangle
        dist_point_interface(current_position, interface, alphas, triangle);
        current_position = alphas(0)*triangle(0)+alphas(1)*triangle(1)+alphas(2)*triangle(2);
        for ( unsigned k = 0; k < 3; ++k) {
            output(i,k) = current_position(k);
        }
    }

    sensors.getPositions() = output;
    sensors.save(output_filename);

    return 0;
}
