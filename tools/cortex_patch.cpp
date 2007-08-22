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
