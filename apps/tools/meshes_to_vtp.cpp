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

inline std::string
param(const char* base,const unsigned n) {
    std::stringstream res;
    res << base << n;
    return res.str();
}

int main( int argc, char **argv) {

    print_version(argv[0]);

    command_usage("Convert meshes OR geometry into a single VTK/VTP file.");
    const char* geom_filename = command_option("-g",nullptr,"Input geometry");

    const char * input[7];
    for (unsigned i=0; i<7; ++i) {
        const std::string& pname = param("-i",i+1);
        input[i] = command_option(pname.c_str(),nullptr,"Input mesh");
    }

    const char * name[7];
    for (unsigned i=0; i<7; ++i) {
        const std::string& pname = param("-n",i+1);
        const std::string& lname = param("",i+1);
        name[i]  = command_option(pname.c_str(),lname.c_str(),"Mesh name");
    }

    const char* output = command_option("-o" ,nullptr,"Output VTP file");

    if (command_option("-h",nullptr,nullptr))
        return 0;

    if ((!input[0] && !geom_filename) || !output) {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }
    
    Geometry geo;

    if (!geom_filename) {
        Meshes meshes;

        // first only read the number of inputs

        for (unsigned i=0; input[i]!=0; ++i) {
            meshes.emplace_back({ input[i], false, name[i] });
            meshes.back().correct_local_orientation();
        }

        geo.import_meshes(meshes);
    } else {
        geo.read(geom_filename);
    }

    geo.info();

    geo.write_vtp(output);

    return 0;

