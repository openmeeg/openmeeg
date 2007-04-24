#include "vecteur.h"
#include "matrice.h"
#include "symmatrice.h"
#include <iostream>
#include <fstream>
#include "mesh3.h"
#include "vect3.h"

using namespace std;
//using namespace CLMatLib;

int main( int argc, char **argv)
{
    if(argc<4)
    {
        cout << " usage: "<< argv[0] << " input_mesh patches (txt) output_texture (txt)" << endl;
        exit(1);
    }

    // Loading mesh for distributed sources
    mesh input_mesh;
    int ls=(int)strlen(argv[1]);
    if (!strcmp(argv[1]+ls-3,"vtk"))      input_mesh.load_vtk(argv[1]);
    else if (!strcmp(argv[1]+ls-3,"geo")) input_mesh.load_3d(argv[1]);
    else if (!strcmp(argv[1]+ls-3,"tri")) input_mesh.load_tri(argv[1]);

    int nb_points = input_mesh.nbr_pts();
    cout << "Nb of vertices : " << nb_points << endl;

    cout << "Reading patches from  : " << argv[2] << endl;
    matrice patchs(argv[2],'t');

    if(patchs.ncol() != 4)
    {
        cerr << "Nb cols should be 4 not " << patchs.ncol() << endl;
        exit(1);
    }

    cout << "Nb patches : " << patchs.nlin() << endl;

    vect3 *center = new vect3();
    float radius;

    matrice mask(nb_points,patchs.nlin());

    for( unsigned int k = 0; k < patchs.nlin(); k += 1 )
    {
        center->_x() = patchs(k,0);
        center->_y() = patchs(k,1);
        center->_z() = patchs(k,2);
        cout << "-----" << endl;
        cout << "Center : " << *center << endl;
        radius = patchs(k,3);
        cout << "Radius : " << radius << " mm" << endl;

        for( int i = 0; i < nb_points; i += 1 )
        {
            if((input_mesh.vctr(i) - *center).norme() > radius)
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
