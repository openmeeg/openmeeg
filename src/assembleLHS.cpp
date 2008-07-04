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

#if WIN32
#define _USE_MATH_DEFINES
#endif

#include <math.h>

#include "matrice.h"
#include "symmatrice.h"
#include "geometry.h"
#include "operators.h"
#include "assemble.h"

template<class T>
void deflat(T &M, int start, int end, double coef)
{// deflate the matrix
    for(int i=start;i<=end;i++)
    {
        #ifdef USE_OMP
        #pragma omp parallel for
        #endif
        for(int j=i;j<=end;j++)
        {
            M(i,j)+=coef;
        }
    }
}

void assemble_LHS(const Geometry &geo,symmatrice &mat,const int GaussOrder)
{
    int offset=0;

    mat = symmatrice(geo.size());

    for(int c=0;c<geo.nb()-1;c++)
    {
        int offset0=offset;
        int offset1=offset+geo.getM(c).nbPts();
        int offset2=offset+geo.getM(c).nbPts()+geo.getM(c).nbTrgs();
        int offset3=offset+geo.getM(c).nbPts()+geo.getM(c).nbTrgs()+geo.getM(c+1).nbPts();

        //Computing S block first because it's needed for the corresponding N block
        if(c==0) operatorS(geo,c,c,GaussOrder,mat,offset1,offset1);
        operatorS(geo,c+1,c,GaussOrder,mat,offset3,offset1);
        operatorS(geo,c+1,c+1,GaussOrder,mat,offset3,offset3);

        //Computing N block
        if(c==0) operatorN(geo,c,c,GaussOrder,mat,offset0,offset0,offset1,offset1);
        operatorN(geo,c+1,c,GaussOrder,mat,offset2,offset0,offset3,offset1);
        operatorN(geo,c+1,c+1,GaussOrder,mat,offset2,offset2,offset3,offset3);

        //Computing D block
        if(c==0) operatorD(geo,c,c,GaussOrder,mat,offset1,offset0);
        if(c!=geo.nb()-2) operatorD(geo,c+1,c,GaussOrder,mat,offset3,offset0);
        operatorD(geo,c,c+1,GaussOrder,mat,offset1,offset2);
        if(c!=geo.nb()-2) operatorD(geo,c+1,c+1,GaussOrder,mat,offset3,offset2);

        offset=offset2;
    }

    //Block multiplications
    //Because only half the matrix is stored, only the lower part of the matrix is treated
    offset=0;
    double K = 1 / (4*M_PI);
    for(int c=0;c<geo.nb()-1;c++)
    {
        int offset0=offset;
        int offset1=offset+geo.getM(c).nbPts();
        int offset2=offset+geo.getM(c).nbPts()+geo.getM(c).nbTrgs();
        int offset3=offset+geo.getM(c).nbPts()+geo.getM(c).nbTrgs()+geo.getM(c+1).nbPts();
        int offset4=offset+geo.getM(c).nbPts()+geo.getM(c).nbTrgs()+geo.getM(c+1).nbPts()+geo.getM(c+1).nbTrgs();

        //Each operator is scaled with the appropriate constant

        //Column 1
        if(c==0) mult(mat,offset0,offset0,offset1,offset1,(geo.sigma_in(c)+geo.sigma_out(c))*K);
        if(c==0) mult(mat,offset1,offset0,offset2,offset1,-2.0*K);
        mult(mat,offset2,offset0,offset3,offset1,(-geo.sigma_out(c))*K);
        mult(mat,offset3,offset0,offset4,offset1,K);

        //Column 2
        if(c==0) mult(mat,offset1,offset1,offset2,offset2,(1.0/geo.sigma_in(c)+1.0/geo.sigma_out(c))*K);
        mult(mat,offset2,offset1,offset3,offset2,K);
        mult(mat,offset3,offset1,offset4,offset2,(-1.0/geo.sigma_out(c))*K);

        //Column 3
        mult(mat,offset2,offset2,offset3,offset3,(geo.sigma_in(c+1)+geo.sigma_out(c+1))*K);
        mult(mat,offset3,offset2,offset4,offset3,-2.0*K);

        //Column 4
        mult(mat,offset3,offset3,offset4,offset4,(1.0/geo.sigma_in(c+1)+1.0/geo.sigma_out(c+1))*K);

        offset=offset2;
    }
    // Deflate the last diagonal block of new 'mat' :
    int newsize = geo.size()-(geo.getM(geo.nb()-1)).nbTrgs();
    offset = newsize-(geo.getM(geo.nb()-1)).nbPts();
    deflat(mat,offset,newsize-1,mat(offset,offset)/(newsize-offset));

    mat = mat.submat(0,newsize-1);
}

void assemble_SurfToVol(const Geometry &geo,matrice &mat,const matrice &points)
{
 // only consider innermost surface points and triangles
  // (for the moment SurfToVol only works for the innermost surface and volume)
  int c=0;
  int offset=0;
  int offset0=offset;
  int offset1=offset+geo.getM(c).nbPts();
  int nbpoints=geo.getM(0).nbPts();
  int nbtriangles = geo.getM(0).nbTrgs();
  double K = 1/(4*M_PI);
  std::cout<<" nbpoints= " << nbpoints <<std::endl;
  std::cout<<" nbtriangles= " << nbtriangles <<std::endl;
  std::cout<< "observation points: " << points.nlin() << std::endl;
   mat = matrice(points.nlin(),nbpoints+nbtriangles);
     // compute S blocks
       operatorSinternal(geo,c,mat,offset1,points);
       mult2(mat,offset0,offset1,offset0+points.nlin(),offset1+geo.getM(0).nbTrgs(),K);
      // compute D blocks
      operatorDinternal(geo,c,mat,offset0,points);
      mult2(mat,offset0,offset0,offset0+points.nlin(),offset1,-(1.0/geo.sigma_in(0))*K);
}
LHS_matrice::LHS_matrice (const Geometry &geo, const int GaussOrder) {
  assemble_LHS(geo,*this,GaussOrder);
}

SurfToVol_matrice::SurfToVol_matrice (const Geometry &geo, const matrice &points) {
  assemble_SurfToVol(geo,*this,points);
}
