/*
Project Name : OpenMEEG

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

namespace OpenMEEG {

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
        mult(mat,nVertexFirstLayer,0,nVertexFirstLayer+nFacesFirstLayer-1,nVertexSources-1,(-1.0/s1i)*K);
        mult(mat,0,0,nVertexFirstLayer-1,nVertexSources-1,K);
    }

    SurfSourceMat::SurfSourceMat (const Geometry &geo, const Mesh& sources, const int GaussOrder) {
        assemble_SurfSourceMat(*this,geo,sources,GaussOrder);
    }

    void assemble_DipSourceMat(Matrix &rhs,const Geometry &geo,vector<Vect3> Rs,vector<Vect3> Qs,const int GaussOrder, const bool adapt_rhs)
    {
        int newsize=geo.size()-(geo.getM(geo.nb()-1)).nbTrgs();
        rhs = Matrix(newsize, Qs.size());
        int nVertexFirstLayer=geo.getM(0).nbPts();

        int nFacesFirstLayer=geo.getM(0).nbTrgs();
        double K = 1.0 / (4*M_PI);

        // First block is nVertexFistLayer
        rhs.set(0);
        Vector prov(rhs.nlin());
        for (size_t s=0; s<Qs.size(); s++) {
            PROGRESSBAR(s,Qs.size());
            prov.set(0);
            operatorDipolePotDer(Rs[s],Qs[s],geo.getM(0),prov,0,GaussOrder,adapt_rhs);
            // Second block is nFaceFistLayer
            if(geo.nb()>1) {
                operatorDipolePot(Rs[s],Qs[s],geo.getM(0),prov,nVertexFirstLayer,GaussOrder,adapt_rhs);
            }
            for (size_t i=0; i<rhs.nlin(); i++) {              
                rhs(i,s) = prov(i);
            }
        }
        // Blocks multiplication
        double s1i=geo.sigma_in(0);
        for( int i=0; i<nVertexFirstLayer; i++ ) {
            for (unsigned j=0; j<rhs.ncol(); j++) {
                rhs(i,j) *= K;
            }
        }
        if (geo.nb()>1) {
            for( int i=0; i<nFacesFirstLayer; i++ ) {
                for (unsigned j=0; j<rhs.ncol(); j++) {
                    rhs(i+nVertexFirstLayer, j) *= (-K/s1i);
                }
            }
        }
    }

    DipSourceMat::DipSourceMat (const Geometry &geo, const Matrix& dipoles, const int GaussOrder, const bool adapt_rhs=true) {
        vector<Vect3> Rs, Qs;

        // Assembling Matrix from discretization :
        unsigned int nd = (unsigned int) dipoles.nlin();
        for( unsigned int i=0; i<nd; i++ )
        {
            Vect3 r(3),q(3);
            for(int j=0;j<3;j++) r(j)   = dipoles(i,j);
            for(int j=3;j<6;j++) q(j-3) = dipoles(i,j);
            Rs.push_back(r); Qs.push_back(q);
        }

        assemble_DipSourceMat(*this,geo,Rs,Qs,GaussOrder,adapt_rhs);
    }

    void assemble_EITSourceMat(Matrix &mat, const Geometry &geo, Matrix& triangleArea, const int GaussOrder)
    {
        // a Matrix to be applied to the scalp-injected current (modulo multiplicative constants)
        // to obtain the Source Term of the EIT foward problem
        int totalsize=geo.size();
        int sourcesize = (geo.getM(geo.nb()-1)).nbTrgs();
        int newsize=totalsize-sourcesize;
        mat=Matrix(newsize,sourcesize);
        // transmat = a big  Matrix of which mat = part of its transpose
        SymMatrix transmat(newsize+sourcesize);
        // airemat = a Matrix to store the surface of triangles on the scalp, for normalizing the injected current
        SymMatrix transairescalp(newsize+sourcesize);
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
        operatorS(geo.getM(c+1),geo.getM(c),transmat,offset3,offset1,GaussOrder);
        mult(transmat,offset3,offset1,offset4,offset2,K/geo.sigma_in(c+1));
        // first compute D, then it will be transposed
        operatorD(geo.getM(c+1),geo.getM(c),transmat,offset3,offset0,GaussOrder);
        mult(transmat,offset3,offset0,offset4,offset1,-K);
        operatorD(geo.getM(c+1),geo.getM(c+1),transmat,offset3,offset2,GaussOrder);
        mult(transmat,offset3,offset2,offset4,offset3,-2.0*K);
        operatorP1P0(geo.getM(c+1), transmat,offset3,offset2);
        mult(transmat,offset3,offset2,offset4,offset3,-1/2.0);
        operatorP1P0(geo.getM(c+1), transairescalp,offset3,offset2);
        // extracting the transpose of the last block of lines of transmat
        // transposing the Matrix
        for(int i=0;i<newsize;i++) {
            for(int j=0;j<sourcesize;j++) {
                mat(i,j) = transmat(newsize+j,i);
                triangleArea(i,j) = transairescalp(newsize+j,i);
            }
        }
    }

    EITSourceMat::EITSourceMat (const Geometry &geo, Matrix& triangleArea, const int GaussOrder) {
        assemble_EITSourceMat(*this,geo, triangleArea, GaussOrder);
    }

}
