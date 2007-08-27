#define _USE_MATH_DEFINES
#include <math.h>

#include "operateurs.h"
#define MU0 1 //1.25e-6


// geo = geometry 
// mat = storage for ferguson matrix
// pts = where the magnetic field is to be computed 
// n   = numbers of places where magnetic field is to be computed
void assemble_ferguson(const Geometry &geo,matrice &mat, const Vect3 *pts,const int n)
{
    int offsetJ=0;
    // Computation of blocks of Ferguson's matrix
    for(int c=0;c<geo.nb();c++)
    {
        int offsetI=0;
        for (int p=0;p<n;p++)
        {
            operateurFerguson(pts[p],geo.getM(c),mat,offsetI,offsetJ);
            offsetI+=3;
        }
        offsetJ+=geo.getM(c).nbPts();
    }

    //Blocks multiplications
    offsetJ=0;
    for(int c=0;c<geo.nb();c++)
    {
        //int offsetI=0;
        mult2(mat,0,offsetJ,mat.nlin(),offsetJ+geo.getM(c).nbPts(),(geo.sigma_in(c)-geo.sigma_out(c))*MU0/(4*M_PI));
        offsetJ+=geo.getM(c).nbPts();
    }

}


void compute_Binf ( const Mesh& sources_mesh, const matrice& squids_positions, matrice& field_at_squids)
// squids positions are layed line by line in the positions matrix (though it is nSquidsX3 matrix)
// field_at_squids is a matrice such that field_at_squids*sources_intensity=Binf
{
    int nSquids=(int)squids_positions.nlin();
    //int nVertices=(int)sources_mesh.nbPts();

    for(int iSquid=0;iSquid<nSquids;iSquid++)
    {
        Vect3 x( squids_positions(iSquid,0), squids_positions(iSquid,1), squids_positions(iSquid,2));
        operateurFerguson(x,sources_mesh,field_at_squids,3*iSquid,0);
    }
}
