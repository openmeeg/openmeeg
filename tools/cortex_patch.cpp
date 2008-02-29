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

#include <iostream>
#include <fstream>
#include <cstring>

#include "vecteur.h"
#include "matrice.h"
#include "symmatrice.h"
#include "mesh3.h"
#include "vect3.h"

using namespace std;

int main( int argc, char **argv)
{
    if(argc<4)
    {
        cout << " usage: "<< argv[0] << " input_mesh patches (txt) output_texture (txt)" << endl;
        exit(1);
    }

    // Loading mesh for distributed sources
    Mesh input_mesh;
    input_mesh.load(argv[1]);

    int nb_points = input_mesh.nbPts();
    cout << "Nb of vertices : " << nb_points << endl;

    cout << "Reading patches from  : " << argv[2] << endl;
    matrice patchs(argv[2],'t');

    if(patchs.ncol() != 4)
    {
        cerr << "Nb cols should be 4 not " << patchs.ncol() << endl;
        exit(1);
    }

    cout << "Nb patches : " << patchs.nlin() << endl;

    Vect3 *center = new Vect3();
    float radius;

    matrice mask(nb_points,patchs.nlin());

    for( unsigned int k = 0; k < patchs.nlin(); k += 1 )
    {
        center[0] = patchs(k,0);
        center[1] = patchs(k,1);
        center[2] = patchs(k,2);
        cout << "-----" << endl;
        cout << "Center : " << *center << endl;
        radius = float(patchs(k,3));
        cout << "Radius : " << radius << " mm" << endl;

        for( int i = 0; i < nb_points; i += 1 )
        {
            if((input_mesh.getPt(i) - *center).norme() > radius)
                mask(i,k) = 0;
            else
                mask(i,k) = 1;
        }
    }

    mask.saveTxt(argv[3]);
    cout << "Mask written in : " << argv[3] << endl;

    delete center;

    return 0;
}
