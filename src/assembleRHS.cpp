#include <vector>
#if WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#include "vecteur.h"
#include "matrice.h"
#include "operateurs.h"

using namespace std;

void assemble_RHS(const Geometry &geo,const Mesh& sources,matrice &mat,const int GaussOrder)
{
    unsigned nVertexSources=sources.nbPts();
    unsigned nVertexFirstLayer=geo.getM(0).nbPts();
    unsigned nFacesFirstLayer=geo.getM(0).nbTrgs();
    cout << endl << "assemble RHS with " << nVertexFirstLayer << " vertices (Sources)" << endl << endl;

    // First block is nVertexFistLayer*nVertexSources
    operateurN(geo.getM(0),sources,mat,0,0,GaussOrder);
    
    // Second block is nFacesFistLayer*nVertexSources
    operateurD(geo.getM(0),sources,mat,(int)nVertexFirstLayer,0,GaussOrder);

    double K=1.0/(4.0*M_PI);

    // First block*=(-1/sigma_inside)
    double s1i=geo.sigma_in(0);
    mult2(mat,nVertexFirstLayer,0,nVertexFirstLayer+nFacesFirstLayer-1,nVertexSources-1,(-1.0/s1i)*K );
    mult2(mat,0,0,nVertexFirstLayer-1,nVertexSources-1,K );
}

void assemble_RHS2(const Geometry &geo,const Mesh &sources, matrice &mat,const int GaussOrder)
{
    unsigned nVertexSources=sources.nbPts();
    unsigned nVertexFirstLayer=geo.getM(0).nbPts();
    unsigned nFacesFirstLayer=geo.getM(0).nbTrgs();

    // S block (nFacesfirstLayer*nFacesSources) is computed in order to speed-up the computation of the N block
    operateurS(geo.getM(0),sources,mat,(int)nVertexFirstLayer,(int)nVertexSources,GaussOrder);

    // First block is nVertexFistLayer*nVertexSources
    operateurN(geo.getM(0),sources,mat,0,0,GaussOrder,(int)nVertexFirstLayer,(int)nVertexSources);

    // Second block is nFacesFistLayer*nVertexSources
    operateurD(geo.getM(0),sources,mat,(int)nVertexFirstLayer,0,GaussOrder);

    double K=1.0/(4.0*M_PI);

    // First block*=(-1/sigma_inside)
    double s1i=geo.sigma_in(0);
    mult2(mat,nVertexFirstLayer,0,nVertexFirstLayer+nFacesFirstLayer-1,nVertexSources-1,(-1.0/s1i)*K );
    mult2(mat,0,0,nVertexFirstLayer-1,nVertexSources-1,K );
}

void assemble_RHS_dipoles(const Geometry &geo,vector<Vect3> Rs,vector<Vect3> Qs,matrice &rhs,const int GaussOrder)
{

    unsigned nVertexFirstLayer=geo.getM(0).nbPts();
    unsigned nFacesFirstLayer=geo.getM(0).nbTrgs();

    double K=1.0/(4*M_PI);

    // First block is nVertexFistLayer
    rhs.set(0);
    for( unsigned s=0; s<Qs.size(); s++ )
    {
        vecteur prov(rhs.nlin());

        prov.set(0);

        operateurDipolePotDer(Rs[s],Qs[s],geo.getM(0),prov,0,GaussOrder);

        // Second block is nFaceFistLayer
        operateurDipolePot(Rs[s],Qs[s],geo.getM(0),prov,nVertexFirstLayer,GaussOrder);

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

// Gradient
void assemble_RHS_dipoles_grad(const Geometry &geo,vector<Vect3> Rs,vector<Vect3> Qs,matrice &rhs,const int GaussOrder)
{

	unsigned int nd=Qs.size();

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

        operateurDipolePotDerGrad(Rs[s],Qs[s],geo.getM(0),prov,0,GaussOrder);

        // Second block is nFaceFistLayer
        operateurDipolePotGrad(Rs[s],Qs[s],geo.getM(0),prov,nVertexFirstLayer,GaussOrder);

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

