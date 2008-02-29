/* FILE: $Id$ */

/*
Project Name : OpenMEEG

author            : $Author$
version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

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
    command_usage("Convert matrices between different formats");
    const char *input_filename = command_option("-i",(const char *) NULL,"Input Matrice");
    const char *output_filename = command_option("-o",(const char *) NULL,"Output Matrice");
    const char *use_symmetric = command_option("-sym",(const char *) NULL,"Matrices are symmetric");
    const char *use_sparse = command_option("-sparse",(const char *) NULL,"Matrices are sparse");
    const char *use_binary = command_option("-bin",(const char *) NULL,"Input matrice is in binary format");
    const char *use_txt = command_option("-txt",(const char *) NULL,"Input matrices is in ascii format");
    const char *use_mat = command_option("-mat",(const char *) NULL,"Input matrices is in matlab format");
    if (command_option("-h",(const char *)0,0)) return 0;

    if (use_symmetric) {
        symmatrice M;
        if (use_binary) {
            M.loadBin(input_filename);
        } else if (use_txt) {
            M.loadTxt(input_filename);
        } else if (use_mat) {
            std::cerr << "Matlab format not supported for symmetric matrices" << std::endl;
        } else {
            M.load(input_filename);
        }
        M.save(output_filename);
    } else if (use_sparse) {
        sparse_matrice M;
        if (use_binary) {
            M.loadBin(input_filename);
        } else if (use_txt) {
            M.loadTxt(input_filename);
        } else if (use_mat) {
            M.loadMat(input_filename);
        } else {
            M.load(input_filename);
        }
        M.save(output_filename);
    } else {
        matrice M;
        if (use_binary) {
            M.loadBin(input_filename);
        } else if (use_txt) {
            M.loadTxt(input_filename);
        } else if (use_mat) {
            M.loadMat(input_filename);
        } else {
            M.load(input_filename);
        }
        M.save(output_filename);
    }

    return 0;
}
