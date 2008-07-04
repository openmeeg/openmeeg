/* FILE: $Id: matrix_convert.cpp 208 2008-02-29 13:28:33Z gramfort $ */

/*
Project Name : OpenMEEG

author            : $Author: gramfort $
version           : $Revision: 208 $
last revision     : $Date: 2008-02-29 14:28:33 +0100 (Ven, 29 fév 2008) $
modified by       : $LastChangedBy: gramfort $
last modified     : $LastChangedDate: 2008-02-29 14:28:33 +0100 (Ven, 29 fév 2008) $

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
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

#include "symmatrice.h"
#include "matrice.h"
#include "sparse_matrice.h"
#include "fast_sparse_matrice.h"

#include "options.h"

using namespace std;

int main( int argc, char **argv) {
    command_usage("Convert symmetric matrices between different formats");
    const char *input_filename = command_option("-i",(const char *) NULL,"Input full matrice");
    const char *output_filename = command_option("-o",(const char *) NULL,"Output full matrice");
    const char *input_format = command_option("-if",(const char *) NULL,"Input file format : ascii, binary, old_binary (should be avoided)");
    const char *output_format = command_option("-of",(const char *) NULL,"Output file format : ascii, binary, old_binary (should be avoided)");
    if (command_option("-h",(const char *)0,0)) return 0;

    if(argc<2 || !input_filename || !output_filename) {
        cout << "Not enough arguments, try the -h option" << endl;
        return 1;
    }

    symmatrice M;
    Maths::ifstream ifs(input_filename);
    Maths::ofstream ofs(output_filename);

    try
    {
        if(input_format) {
            ifs >> Maths::format(input_format) >> M;
        } else {
            ifs >> M;
        }

        if(output_format) {
            ofs << Maths::format(output_format) << M;
        } else {
            ofs << Maths::format(output_filename,Maths::format::FromSuffix) << M;
        }
    } catch (std::string s) {
        std::cerr << s << std::endl;
    }

    return 0;
}
