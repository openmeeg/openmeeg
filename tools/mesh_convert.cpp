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

#include "mesh3.h"
#include "options.h"

using namespace std;

int main( int argc, char **argv)
{
    command_usage("Convert mesh between different formats");
    const char *input_filename = command_option("-i",(const char *) NULL,"Input Mesh");
    const char *output_filename = command_option("-o",(const char *) NULL,"Output Mesh");
    const double tx = command_option("-tx",0.0,"Translation along the x axis");
    const double ty = command_option("-ty",0.0,"Translation along the y axis");
    const double tz = command_option("-tz",0.0,"Translation along the z axis");
    const double vx = command_option("-vx",1.0,"Scaling along the x axis");
    const double vy = command_option("-vy",1.0,"Scaling along the y axis");
    const double vz = command_option("-vz",1.0,"Scaling along the z axis");
    const bool apply_asa_flip = command_option("-flip",false,"Rotating axis if mesh comes from ASA");
    if (command_option("-h",(const char *)0,0)) return 0;

    Mesh M;
    M.load(input_filename,false);

    for( unsigned int i = 0; i < unsigned(M.nbPts()); i += 1 )
    {
        Vect3& pt = M[i];
        if (apply_asa_flip) {
            double tmp;
            tmp = pt(0);
            pt(0) = pt(1);
            pt(1) = tmp;
            pt(2) = -pt(2);
        }
        pt(0) = pt(0)+tx*vx;
        pt(1) = pt(1)+ty*vy;
        pt(2) = pt(2)+tz*vz;
    }

    M.save(output_filename);

    return 0;
}
