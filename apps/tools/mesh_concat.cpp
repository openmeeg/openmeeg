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

using namespace std;
using namespace OpenMEEG;

int main( int argc, char **argv)
{
    print_version(argv[0]);

    command_usage("Concat 2 mesh and save the result");
    
    const char *input_filename1 = command_option("-i1", (const char *) "", "Input Mesh 1");
    const char *input_filename2 = command_option("-i2", (const char *) "", "Input Mesh 2");
    const char *output_filename = command_option("-o", (const char *) "", "Output Mesh");

    if ( command_option("-h", (const char *)0, 0)) return 0;

    if ( argc < 2 ) {
        cout << "Not enough arguments, try the -h option" << endl;
        return 1;
    }

    Mesh m1(input_filename1);

    Mesh m2(input_filename2);
    Mesh m3;

    m3.merge(m1, m2);
    m3.save(output_filename);

    return 0;
}
