#include <vector>
#if WIN32
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#include "vecteur.h"
#include "matrice.h"
#include "operateurs.h"
#include "analytiques.h"

using namespace std;

void assemble_RHSmatrix( geometry &geo, mesh &sources, matrice &mat)
{
    unsigned nVertexSources=sources.nbr_pts();
    unsigned nVertexFirstLayer=geo.getM(0).nbr_pts();
    unsigned nFacesFirstLayer=geo.getM(0).nbr_trg();
    cout << endl << "assemble RHS with " << nVertexFirstLayer << " vertices (Sources)" << endl << endl;

    // First block is nVertexFistLayer*nVertexSources
    operateurN( geo.getM(0), sources, mat, 0, 0 );

    // Second block is nFacesFistLayer*nVertexSources
    operateurD( geo.getM(0), sources, mat,  (int)nVertexFirstLayer, 0);

    double K=1.0/(4.0*M_PI);

    // First block*=(-1/sigma_inside)
    double s1i=geo.sigma_in(0);
    mult2(mat,nVertexFirstLayer,0,nVertexFirstLayer+nFacesFirstLayer-1,nVertexSources-1,(-1.0/s1i)*K );
    mult2(mat,0,0,nVertexFirstLayer-1,nVertexSources-1,K );
}

void assemble_RHS2matrix( geometry &geo, mesh &sources, matrice &mat)
{
    unsigned nVertexSources=sources.nbr_pts();
    unsigned nVertexFirstLayer=geo.getM(0).nbr_pts();
    unsigned nFacesFirstLayer=geo.getM(0).nbr_trg();

    // S block (nFacesfirstLayer*nFacesSources) is computed in order to speed-up the computation of the N block
    operateurS(geo.getM(0),sources,mat,(int)nVertexFirstLayer, (int)nVertexSources);

    // First block is nVertexFistLayer*nVertexSources
    operateurN( geo.getM(0), sources, mat, 0, 0, (int)nVertexFirstLayer, (int)nVertexSources);

    // Second block is nFacesFistLayer*nVertexSources
    operateurD( geo.getM(0), sources, mat,  (int)nVertexFirstLayer, 0);

    double K=1.0/(4.0*M_PI);

    // First block*=(-1/sigma_inside)
    double s1i=geo.sigma_in(0);
    mult2(mat,nVertexFirstLayer,0,nVertexFirstLayer+nFacesFirstLayer-1,nVertexSources-1,(-1.0/s1i)*K );
    mult2(mat,0,0,nVertexFirstLayer-1,nVertexSources-1,K );
}

void assemble_RHS_dipoles_matrice( geometry &geo, vector<vect3> Rs, vector<vect3> Qs, matrice &rhs)
{

    unsigned nVertexFirstLayer=geo.getM(0).nbr_pts();
    unsigned nFacesFirstLayer=geo.getM(0).nbr_trg();

    double K=1.0/(4*M_PI);

    // First block is nVertexFistLayer
    rhs.set(0);
    for( unsigned s=0; s<Qs.size(); s++ ) 
    {
        vecteur prov(rhs.nlin());

        prov.set(0);

        operateurDipolePotDer(Rs[s],Qs[s],geo.getM(0),prov,0);

        // Second block is nFaceFistLayer
        operateurDipolePot(Rs[s],Qs[s],geo.getM(0),prov,nVertexFirstLayer);

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

void assemble_RHSvector( geometry &geo, vector<vect3> Rs, vector<vect3> Qs, vecteur &rhs)
{
    unsigned nVertexFirstLayer=geo.getM(0).nbr_pts();
    unsigned nFacesFirstLayer=geo.getM(0).nbr_trg();

    double K=1.0/(4*M_PI);

    // First block is nVertexFistLayer
    rhs.set(0);
    for( unsigned s=0; s<Qs.size(); s++ )
    {
        operateurDipolePotDer(Rs[s],Qs[s],geo.getM(0),rhs,0);
    }

    // Second block is nFaceFistLayer
    for( unsigned s=0;s<Qs.size();s++ )
    {
        operateurDipolePot(Rs[s],Qs[s],geo.getM(0),rhs,nVertexFirstLayer);
    }
    
    // Blocks multiplication
    double s1i=geo.sigma_in(0);
    for(unsigned i=0;i<nVertexFirstLayer;i++)
    {
        rhs(i)*=K;
    }
    for(unsigned i=0;i<nFacesFirstLayer;i++)
    {
        rhs(i+nVertexFirstLayer)*=(-K/s1i);
    }
}
