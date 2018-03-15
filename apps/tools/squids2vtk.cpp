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

#include <iostream>
#include "options.h"
#include "matrix.h"
#include "symmatrix.h"
#include "vector.h"
#include "om_utils.h"

using namespace OpenMEEG;

int main( int argc, char** argv)
{
    print_version(argv[0]);

    command_usage("Convert squids in text file to a vtk file for vizualisation");
    const char *input_filename = command_option("-i",(const char *) NULL,"Squids positions in original coordinate system");
    const char *output_filename = command_option("-o",(const char *) NULL,"Squids positions with orientations in vtk format");
    if (command_option("-h",(const char *)0,0)) return 0;

    if(!input_filename || !output_filename) {
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

    for (unsigned int i = 0; i < squids.nlin();++i)
        ofs << squids(i,0) << ' ' << squids(i,1) << ' ' << squids(i,2) << std::endl;

    ofs << "POINT_DATA " << squids.nlin() << std::endl
        << "NORMALS normals float" << std::endl;

    for (unsigned int i = 0; i < squids.nlin();++i)
        ofs << squids(i,3) << ' ' << squids(i,4) << ' ' << squids(i,5) << std::endl;

    return 0;
}
