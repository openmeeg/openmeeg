/* FILE: $Id: sparse_#include "Maths.H" 235 2008-04-30 13:39:08Z gramfort $ */

/*
Project Name : OpenMEEG

author            : $Author: $
version           : $Revision: $
last revision     : $Date: $
modified by       : $LastChangedBy: $
last modified     : $LastChangedDate: $

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

#include "mesh3.h"
#include "options.h"
#include "matrice_dcl.h"

using namespace std;

int main( int argc, char **argv)
{
    command_usage("Convert mesh file to a dipole file");
    const char *input_filename = command_option("-i",(const char *) NULL,"Input Mesh");
    const char *output_filename = command_option("-o",(const char *) NULL,"Output .dip file");
    if (command_option("-h",(const char *)0,0)) return 0;

    if(!input_filename || !output_filename) {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    Mesh M;
    M.load(input_filename,false);

    matrice mat(M.nbPts(),6);
    for(int i = 0; i < M.nbPts(); ++i)
    {
        mat(i,0) = M.getPt(i).x();
        mat(i,1) = M.getPt(i).y();
        mat(i,2) = M.getPt(i).z();
        mat(i,3) = M.normal(i).x();
        mat(i,4) = M.normal(i).y();
        mat(i,5) = M.normal(i).z();
    }
    mat.saveTxt(output_filename);

    return 0;
}
