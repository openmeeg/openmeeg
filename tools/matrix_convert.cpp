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

#include "symmatrix.h"
#include "matrix.h"
#include "sparse_matrix.h"
#include "fast_sparse_matrix.h"

#include "options.h"

using namespace std;
using namespace OpenMEEG;

int main( int argc, char **argv)
{
    print_version(argv[0]);

    command_usage("Convert full matrices between different formats");
    const char *input_filename = command_option("-i",(const char *) NULL,"Input full Matrix");
    const char *output_filename = command_option("-o",(const char *) NULL,"Output full Matrix");
    const char *input_format = command_option("-if",(const char *) NULL,"Input file format : ascii, binary, tex, matlab");
    const char *output_format = command_option("-of",(const char *) NULL,"Output file format : ascii, binary, tex, matlab");
    if (command_option("-h",(const char *)0,0)) return 0;

    if(argc<2 || !input_filename || !output_filename) {
        cout << "Not enough arguments, try the -h option" << endl;
        return 1;
    }

    Matrix M;
    maths::ifstream ifs(input_filename);
    maths::ofstream ofs(output_filename);

    try
    {
        if(input_format) {
            ifs >> maths::format(input_format) >> M;
        } else {
            ifs >> M;
        }

        if(output_format) {
            ofs << maths::format(output_format) << M;
        } else {
            ofs << maths::format(output_filename,maths::format::FromSuffix) << M;
        }
    } catch (std::string s) {
        std::cerr << s << std::endl;
    }

    return 0;
}
