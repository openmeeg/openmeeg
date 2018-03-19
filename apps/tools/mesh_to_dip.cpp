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

#include "mesh.h"
#include "options.h"
#include "matrix.h"

using namespace std;
using namespace OpenMEEG;

int main( int argc, char **argv)
{
    print_version(argv[0]);

    command_usage("Convert mesh file to a dipole file");
    const char *input_filename  = command_option("-i", (const char *) NULL, "Input Mesh");
    const char *output_filename = command_option("-o", (const char *) NULL, "Output .dip file");

    if (command_option("-h",(const char *)0,0)) return 0;

    if(!input_filename || !output_filename) {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    Mesh m(input_filename, false);

    Matrix mat(m.nb_vertices(), 6);

    unsigned i = 0;

    for ( Mesh::const_vertex_iterator vit = m.vertex_begin(); vit != m.vertex_end(); ++vit, ++i)
    {
        mat(i, 0) = (*vit)->x();
        mat(i, 1) = (*vit)->y();
        mat(i, 2) = (*vit)->z();
        Normal n = m.normal(**vit);
        mat(i, 3) = n.x();
        mat(i, 4) = n.y();
        mat(i, 5) = n.z();
    }

    mat.save(output_filename);

    return 0;
}
