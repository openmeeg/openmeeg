/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.sophia.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.sophia.inria.fr)

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
#include <geometry_io.h>
#include <matrix.h>
#include <options.h>

using namespace std;
using namespace OpenMEEG;

int main( int argc, char **argv) {

    print_version(argv[0]);

    command_usage("Convert meshes into a single VTP file.");
    const char * input[6];
    const char * output;
    const char * name[6];
    input[0] = command_option("-i1", (const char *) NULL, "Input Mesh");
    input[1] = command_option("-i2", (const char *) NULL, "Input Mesh");
    input[2] = command_option("-i3", (const char *) NULL, "Input Mesh");
    input[3] = command_option("-i4", (const char *) NULL, "Input Mesh");
    input[4] = command_option("-i5", (const char *) NULL, "Input Mesh");
    input[5] = command_option("-i6", (const char *) NULL, "Input Mesh");
    name[0]  = command_option("-n1", (const char *) NULL, "Mesh Name");
    name[1]  = command_option("-n2", (const char *) NULL, "Mesh Name");
    name[2]  = command_option("-n3", (const char *) NULL, "Mesh Name");
    name[3]  = command_option("-n4", (const char *) NULL, "Mesh Name");
    name[4]  = command_option("-n5", (const char *) NULL, "Mesh Name");
    name[5]  = command_option("-n6", (const char *) NULL, "Mesh Name");
    output   = command_option("-o" , (const char *) NULL, "Output VTP file");

    const char *output_filename = command_option("-o",(const char *) NULL,"Output VTP mesh");

    if ( command_option("-h", (const char *)0, 0) )
        return 0;

    if ( !input[0] || !output_filename ) {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }
    
    Meshes   meshes;

    unsigned nb_vertices = 0;
    unsigned i_input     = 0;

    // first only read the number of inputs
    while ( input[i_input] != 0 ) {
        Mesh m(input[i_input], false, name[i_input]);
        meshes.push_back(m);
        ++i_input;
    }

    Geometry geo;

    geo.import_meshes(meshes);

    geo.info();

    geo.write_vtp(output);

    return 0;
}
