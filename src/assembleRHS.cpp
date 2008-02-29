/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
lastrevision      : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
*/

#include <vector>
#if WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#include "vecteur.h"
#include "matrice.h"
#include "operators.h"
#include "assemble.h"

using namespace std;

void assemble_RHS(matrice &mat,const Geometry &geo,const Mesh& sources,const int GaussOrder)
{
    int newsize = geo.size()-(geo.getM(geo.nb()-1)).nbTrgs();
    mat = matrice(newsize,sources.nbPts());
    mat.set(0.0);

    unsigned nVertexSources=sources.nbPts();
    unsigned nVertexFirstLayer=geo.getM(0).nbPts();
    unsigned nFacesFirstLayer=geo.getM(0).nbTrgs();
    cout << endl << "assemble RHS with " << nVertexSources << " sources" << endl << endl;

    // First block is nVertexFistLayer*nVertexSources
    operatorN(geo.getM(0),sources,mat,0,0,GaussOrder);

    // Second block is nFacesFistLayer*nVertexSources
    operatorD(geo.getM(0),sources,mat,(int)nVertexFirstLayer,0,GaussOrder);

    double K=1.0/(4.0*M_PI);

    // First block*=(-1/sigma_inside)
    double s1i=geo.sigma_in(0);
    mult2(mat,nVertexFirstLayer,0,nVertexFirstLayer+nFacesFirstLayer-1,nVertexSources-1,(-1.0/s1i)*K);
    mult2(mat,0,0,nVertexFirstLayer-1,nVertexSources-1,K);
}

RHS_matrice::RHS_matrice (const Geometry &geo, const Mesh& sources, const int GaussOrder) {
    assemble_RHS(*this,geo,sources,GaussOrder);
}

void assemble_RHSdip(matrice &rhs,const Geometry &geo,vector<Vect3> Rs,vector<Vect3> Qs,const int GaussOrder)
{
    int newsize=geo.size()-(geo.getM(geo.nb()-1)).nbTrgs();
    rhs = matrice(newsize, Qs.size());

    unsigned nVertexFirstLayer=geo.getM(0).nbPts();
    unsigned nFacesFirstLayer=geo.getM(0).nbTrgs();

    double K=1.0/(4*M_PI);

    // First block is nVertexFistLayer
    rhs.set(0);
    for( unsigned s=0; s<Qs.size(); s++ )
    {
        vecteur prov(rhs.nlin());

        prov.set(0);

        operatorDipolePotDer(Rs[s],Qs[s],geo.getM(0),prov,0,GaussOrder);

        // Second block is nFaceFistLayer
        operatorDipolePot(Rs[s],Qs[s],geo.getM(0),prov,nVertexFirstLayer,GaussOrder);

        for (unsigned i=0; i<rhs.nlin(); i++)
        {
            rhs(i, s) = prov(i);
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
    for( unsigned i=0; i<nFacesFirstLayer; i++ )
    {
        for (unsigned j=0; j<rhs.ncol(); j++)
        {
            rhs(i+nVertexFirstLayer, j) *= (-K/s1i);
        }
    }
}

RHSdip_matrice::RHSdip_matrice (const Geometry &geo, vector<Vect3> Rs, vector<Vect3> Qs, const int GaussOrder) {
    assemble_RHSdip(*this,geo,Rs,Qs,GaussOrder);
}

// Gradient
void assemble_RHSdip_grad(matrice &rhs,const Geometry &geo,vector<Vect3> Rs,vector<Vect3> Qs,const int GaussOrder)
{
    unsigned int nd=Qs.size();

    int newsize = geo.size()-(geo.getM(geo.nb()-1)).nbTrgs();
    rhs = matrice(newsize, 6*nd); // 6 derivatives!

    unsigned nVertexFirstLayer=geo.getM(0).nbPts();
    unsigned nFacesFirstLayer=geo.getM(0).nbTrgs();

    double K=1.0/(4*M_PI);

    // First block is nVertexFistLayer
    rhs.set(0);
    for( unsigned s=0; s<nd; s++ )
    {
        vecteur prov[6];
        for (int d=0;d<6;d++) {
            prov[d]=vecteur(rhs.nlin());
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

RHSdip_grad_matrice::RHSdip_grad_matrice (const Geometry &geo, vector<Vect3> Rs, vector<Vect3> Qs, const int GaussOrder) {
    assemble_RHSdip_grad(*this,geo,Rs,Qs,GaussOrder);
}

void assemble_EITsource(const Geometry &geo, matrice &mat, matrice &airescalp, const int GaussOrder)
{
// Une matrice qu'il suffira de multiplier par la valeur du courant
// injecte sur le scalp (modulo les constantes multiplicatives)
//  pour obtenir le second membre du pb direct de l'EIT
// 
    int newtaille = mat.nlin();
    int sourcetaille = mat.ncol();
// transmat = une grande matrice dont mat = une partie de la transposee
    symmatrice transmat(newtaille+sourcetaille);
// airemat = une matrice qui va servir a stocker l'aire des triangles du 
// scalp, pour l'injection du courant
    symmatrice transairescalp(newtaille+sourcetaille);
    int c;
    int offset=0;
    int offset0;
    int offset1;
    int offset2;
    int offset3;
    int offset4;
    double K=1.0/(4*M_PI);

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

// calcul de S 
    operatorS(geo,c,c-1,GaussOrder,transmat,offset3,offset1);
    mult(transmat,offset3,offset1,offset4,offset2,-1.0*K);
// on calcule d'abord D puis on transposera
    operatorD(geo,c+1,c,GaussOrder,transmat,offset3,offset0);
    mult(transmat,offset3,offset0,offset4,offset1,K);
    operatorD(geo,c+1,c+1,GaussOrder,transmat,offset3,offset2);
    mult(transmat,offset3,offset2,offset4,offset3,-1.0*K);
    operatorP1P0(geo,c+1, transmat,offset3,offset2);
    operatorP1P0(geo,c+1, transairescalp,offset3,offset2);
// on extrait la transposee du dernier bloc de lignes de transmat
    std::cout<<"offset0 "<<offset0<<std::endl;
    std::cout<<"offset1 "<<offset1<<std::endl;
    std::cout<<"offset2 "<<offset2<<std::endl;
    std::cout<<"offset3 "<<offset3<<std::endl;

// on transpose la matrice
    std::cout<<"last element "<<transmat(newtaille+sourcetaille-1,newtaille-1)<<std::endl;
        for(int i=0;i<newtaille;i++) 
            for(int j=0;j<sourcetaille;j++) {
                mat(i,j) = transmat(newtaille+j,i);
                airescalp(i,j) = 2*transairescalp(newtaille+j,i);
            }
}

