/* FILE: $Id: assembleRHS.cpp 257 2008-07-08 16:03:45Z gramfort $ */

/*
Project Name : OpenMEEG

version           : $Revision: 257 $
last revision     : $Date: 2008-07-08 18:03:45 +0200 (Tue, 08 Jul 2008) $
modified by       : $LastChangedBy: gramfort $
last modified     : $LastChangedDate: 2008-07-08 18:03:45 +0200 (Tue, 08 Jul 2008) $

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

#include <vector>
#if WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#include "vector.h"
#include "matrix.h"
#include "operators.h"
#include "assemble.h"
#include <fstream>
using namespace std;

void assemble_SurfSourceMat(Matrix &mat,const Geometry &geo,const Mesh& sources,const int GaussOrder)
{
    int newsize = geo.size()-(geo.getM(geo.nb()-1)).nbTrgs();
    mat = Matrix(newsize,sources.nbPts());
    mat.set(0.0);

    unsigned nVertexSources=sources.nbPts();
    unsigned nVertexFirstLayer=geo.getM(0).nbPts();
    unsigned nFacesFirstLayer=geo.getM(0).nbTrgs();
    cout << endl << "assemble SurfSourceMat with " << nVertexSources << " sources" << endl << endl;

    // First block is nVertexFistLayer*nVertexSources
    operatorN(geo.getM(0),sources,mat,0,0,GaussOrder);

    // Second block is nFacesFistLayer*nVertexSources
    operatorD(geo.getM(0),sources,mat,(int)nVertexFirstLayer,0,GaussOrder);

    double K = 1.0 / (4.0*M_PI);

    // First block*=(-1/sigma_inside)
    double s1i=geo.sigma_in(0);
    mult2(mat,nVertexFirstLayer,0,nVertexFirstLayer+nFacesFirstLayer-1,nVertexSources-1,(-1.0/s1i)*K);
    mult2(mat,0,0,nVertexFirstLayer-1,nVertexSources-1,K);
}

SurfSourceMat::SurfSourceMat (const Geometry &geo, const Mesh& sources, const int GaussOrder) {
    assemble_SurfSourceMat(*this,geo,sources,GaussOrder);
}

void assemble_DipSourceMat(Matrix &rhs,const Geometry &geo,vector<Vect3> Rs,vector<Vect3> Qs,const int GaussOrder)
{ 
    int newsize=geo.size()-(geo.getM(geo.nb()-1)).nbTrgs();
    rhs = Matrix(newsize, Qs.size());
     unsigned nVertexFirstLayer=geo.getM(0).nbPts();

    unsigned nFacesFirstLayer=geo.getM(0).nbTrgs();
    double K = 1.0 / (4*M_PI);

    // First block is nVertexFistLayer
    rhs.set(0);
    for (unsigned s=0; s<Qs.size(); s++)
    {  
      Vector prov(rhs.nlin());
      prov.set(0);
        operatorDipolePotDer(Rs[s],Qs[s],geo.getM(0),prov,0,GaussOrder);
        // Second block is nFaceFistLayer
	if(geo.nb()>1){
        operatorDipolePot(Rs[s],Qs[s],geo.getM(0),prov,nVertexFirstLayer,GaussOrder);
	
	}
        for (unsigned i=0; i<rhs.nlin(); i++)
	  {              rhs(i,s) = prov(i);
        }
    }
    // Blocks multiplication
    double s1i=geo.sigma_in(0);
    for( unsigned i=0; i<nVertexFirstLayer; i++ )
    {
        for (unsigned j=0; j<rhs.ncol(); j++) {
            rhs(i,j) *= K;
        }
    }
    if(geo.nb()>1){
    for( unsigned i=0; i<nFacesFirstLayer; i++ )
    {
        for (unsigned j=0; j<rhs.ncol(); j++)
	  { 
            rhs(i+nVertexFirstLayer, j) *= (-K/s1i);
        }
    }
    }
}

DipSourceMat::DipSourceMat (const Geometry &geo, vector<Vect3> Rs, vector<Vect3> Qs, const int GaussOrder) {
    assemble_DipSourceMat(*this,geo,Rs,Qs,GaussOrder);
}

// Gradient
void assemble_DipSourceGradMat(Matrix &rhs,const Geometry &geo,vector<Vect3> Rs,vector<Vect3> Qs,const int GaussOrder)
{
    unsigned int nd=Qs.size();

    int newsize = geo.size()-(geo.getM(geo.nb()-1)).nbTrgs();
    rhs = Matrix(newsize, 6*nd); // 6 derivatives!

    unsigned nVertexFirstLayer=geo.getM(0).nbPts();
    unsigned nFacesFirstLayer=geo.getM(0).nbTrgs();

    double K = 1.0 / (4*M_PI);

    // First block is nVertexFistLayer
    rhs.set(0);
    for( unsigned s=0; s<nd; s++ )
    {
        Vector prov[6];
        for (int d=0;d<6;d++) {
            prov[d]=Vector(rhs.nlin());
            prov[d].set(0);
        }

        operatorDipolePotDerGrad(Rs[s],Qs[s],geo.getM(0),prov,0,GaussOrder);

        // Second block is nFaceFistLayer
        operatorDipolePotGrad(Rs[s],Qs[s],geo.getM(0),prov,nVertexFirstLayer,GaussOrder);

        for (unsigned i=0; i<rhs.nlin(); i++)
            for (int d=0;d<6;d++)
                rhs(i, 6*s+d) = prov[d](i);
    }

    // Blocks multiplication
    double s1i=geo.sigma_in(0);
    for( unsigned i=0; i<nVertexFirstLayer; i++ )
    {
        for (unsigned j=0; j<6*nd; j++) {
            rhs(i,j) *= K;
        }
    }
    for( unsigned i=0; i<nFacesFirstLayer; i++ )
    {
        for (unsigned j=0; j<6*nd; j++)
        {
            rhs(i+nVertexFirstLayer, j) *= (-K/s1i);
        }
    }
}

DipSourceGradMat::DipSourceGradMat (const Geometry &geo, vector<Vect3> Rs, vector<Vect3> Qs, const int GaussOrder) {
    assemble_DipSourceGradMat(*this,geo,Rs,Qs,GaussOrder);
}

void assemble_EITsource(const Geometry &geo, Matrix &mat, Matrix &airescalp, const int GaussOrder)
{
// a Matrix to be applied to the scalp-injected current (modulo multiplicative constants)
// to obtain the Source Term of the EIT foward problem 
    int newtaille = mat.nlin();
    int sourcetaille = mat.ncol();
// transmat = a big  Matrix of which mat = part of its transpose
    SymMatrix transmat(newtaille+sourcetaille);
// airemat = a Matrix to store the surface of triangles on the scalp, for normalizing the injected current
    SymMatrix transairescalp(newtaille+sourcetaille);
    int c;
    int offset=0;
    int offset0;
    int offset1;
    int offset2;
    int offset3;
    int offset4;

    double K = 1.0 / (4*M_PI);

   for(c=0;c<geo.nb()-1;c++)
        {
            offset0=offset;
            offset1=offset+geo.getM(c).nbPts();
            offset2=offset+geo.getM(c).nbPts()+geo.getM(c).nbTrgs();
            offset3=offset+geo.getM(c).nbPts()+geo.getM(c).nbTrgs()+geo.getM(c+1).nbPts();
	    offset4=offset+geo.getM(c).nbPts()+geo.getM(c).nbTrgs()+geo.getM(c+1).nbPts()+geo.getM(c+1).nbTrgs();
            offset=offset2;
        }
    c=geo.nb()-2;

// compute S 
    operatorS(geo,c+1,c,GaussOrder,transmat,offset3,offset1);
    mult(transmat,offset3,offset1,offset4,offset2,-1.0*K);
// first compute D, then it will be transposed
    operatorD(geo,c+1,c,GaussOrder,transmat,offset3,offset0);
    mult(transmat,offset3,offset0,offset4,offset1,K);
    operatorD(geo,c+1,c+1,GaussOrder,transmat,offset3,offset2);
    mult(transmat,offset3,offset2,offset4,offset3,-1.0*K);
    operatorP1P0(geo,c+1, transmat,offset3,offset2);
    operatorP1P0(geo,c+1, transairescalp,offset3,offset2);
// extracting the transpose of the last block of lines of transmat
    std::cout<<"offset0 "<<offset0<<std::endl;
    std::cout<<"offset1 "<<offset1<<std::endl;
    std::cout<<"offset2 "<<offset2<<std::endl;
    std::cout<<"offset3 "<<offset3<<std::endl;

// transposing the Matrix
    std::cout<<"last element "<<transmat(newtaille+sourcetaille-1,newtaille-1)<<std::endl;
        for(int i=0;i<newtaille;i++) 
            for(int j=0;j<sourcetaille;j++) {
                mat(i,j) = transmat(newtaille+j,i);
                airescalp(i,j) = 2*transairescalp(newtaille+j,i);
            }
}




