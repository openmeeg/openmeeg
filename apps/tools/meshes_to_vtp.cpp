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

#include <geometry.h>
#include <matrix.h>
#include <options.h>

using namespace std;
using namespace OpenMEEG;

int main( int argc, char **argv) {

    print_version(argv[0]);

    command_usage("Convert meshes OR geometry into a single VTK/VTP file.");
    const char * geom_filename;
    const char * input[7];
    const char * name[7];
    const char * output;
    geom_filename = command_option("-g", (const char *) NULL, "Input geometry");
    input[0] = command_option("-i1", (const char *) NULL, "Input mesh");
    input[1] = command_option("-i2", (const char *) NULL, "Input mesh");
    input[2] = command_option("-i3", (const char *) NULL, "Input mesh");
    input[3] = command_option("-i4", (const char *) NULL, "Input mesh");
    input[4] = command_option("-i5", (const char *) NULL, "Input mesh");
    input[5] = command_option("-i6", (const char *) NULL, "Input mesh");
    input[6] = command_option("-i7", (const char *) NULL, "Input mesh");
    name[0]  = command_option("-n1", (const char *) "1", "Mesh name");
    name[1]  = command_option("-n2", (const char *) "2", "Mesh name");
    name[2]  = command_option("-n3", (const char *) "3", "Mesh name");
    name[3]  = command_option("-n4", (const char *) "4", "Mesh name");
    name[4]  = command_option("-n5", (const char *) "5", "Mesh name");
    name[5]  = command_option("-n6", (const char *) "6", "Mesh name");
    name[6]  = command_option("-n7", (const char *) "7", "Mesh name");
    output   = command_option("-o" , (const char *) NULL, "Output VTP file");

    if ( command_option("-h", (const char *)0, 0) )
        return 0;

    if ( (!input[0] && !geom_filename) || !output ) {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }
    
    Geometry geo;

    if (!geom_filename)
    {
        Meshes meshes;

        unsigned i_input     = 0;

        // first only read the number of inputs
        while ( input[i_input] != 0 ) {
            Mesh m(input[i_input], false, name[i_input]);
            m.correct_local_orientation();
            meshes.push_back(m);
            ++i_input;
        }

        geo.import_meshes(meshes);
    }
    else
    {
        geo.read(geom_filename);
    }

    geo.info();

    geo.write_vtp(output);

    return 0;
}
