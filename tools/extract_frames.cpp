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

#include "matrice.h"
#include "symmatrice.h"
#include "vecteur.h"
#include "options.h"

using namespace std;

int main( int argc, char **argv)
{
    command_usage("Simple tool to select a few time frames in a dataset.");
    const char *input_filename = command_option("-i",(const char *) SRCPATH("data/Computations/Head1/Head1.src"),"Dataset from which frames are extracted");
    const char *output_filename = command_option("-o",(const char *) "extracted_frames.txt","Extracted time frames");
    const size_t first_frame = command_option("-f",0,"Index of first frame");
    const size_t length = command_option("-l",1,"Nb of frames to extract");
    const double mult = command_option("-m",1.0,"Nb of frames to extract");
    if (command_option("-h",(const char *)0,0)) return 0;

    matrice in(input_filename,'t');
    assert((first_frame+length) <= in.ncol());
    matrice out(in.nlin(),length);
    for( unsigned int i = 0; i < length; i += 1 )
    {
        out.setcol(i,in.getcol(first_frame+i));
    }
    if (mult != 1.0) {
        for (size_t j=0;j<out.ncol();j++)
            for (size_t i=0;i<out.nlin();i++)
            {
                out(i,j) = out(i,j)*mult;
            }
    }
    out.saveTxt(output_filename);
    return 0;
}
