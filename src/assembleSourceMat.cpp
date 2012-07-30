/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
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

#include "vector.h"
#include "matrix.h"
#include "danielsson.h"
#include "operators.h"
#include "assemble.h"
#include "sensors.h"
#include <fstream>

namespace OpenMEEG {

    using namespace std;

    void assemble_SurfSourceMat(Matrix& mat,const Geometry& geo,const Mesh& sources,const int gauss_order) {

        const int newsize = geo.size()-geo.getM(geo.nb()-1).nbTrgs();
        mat = Matrix(newsize,sources.nbPts());
        mat.set(0.0);

        const unsigned nVertexSources = sources.nbPts();
        const unsigned nVertexFirstLayer = geo.getM(0).nbPts();
        const unsigned nFacesFirstLayer = geo.getM(0).nbTrgs();
        std::cout << std::endl << "assemble SurfSourceMat with " << nVertexSources << " sources" << std::endl << std::endl;

        //  First block is nVertexFistLayer*nVertexSources.

        operatorN(geo.getM(0),sources,mat,0,0,gauss_order);

        //  Second block is nFacesFistLayer*nVertexSources.

        operatorD(geo.getM(0),sources,mat,static_cast<int>(nVertexFirstLayer),0,gauss_order);

        //  First block*=(-1/sigma_inside).

        const double K = 1.0/(4.0*M_PI);
        const double s1i = geo.sigma_in(0);

        mult(mat,nVertexFirstLayer,0,nVertexFirstLayer+nFacesFirstLayer-1,nVertexSources-1,-K/s1i);
        mult(mat,0,0,nVertexFirstLayer-1,nVertexSources-1,K);
    }

    SurfSourceMat::SurfSourceMat(const Geometry& geo,const Mesh& sources,const int gauss_order) {
        assemble_SurfSourceMat(*this,geo,sources,gauss_order);
    }

    void assemble_DipSourceMat(Matrix& rhs,const Geometry& geo,const Matrix& dipoles,
                               const int gauss_order,const bool adapt_rhs,const bool dipoles_in_cortex)
    {
        const int newsize = geo.size()-(geo.getM(geo.nb()-1)).nbTrgs();
        const size_t n_dipoles = dipoles.nlin();
        rhs = Matrix(newsize,n_dipoles);
        rhs.set(0);

        const double K = 1.0/(4*M_PI);

        //  First block is nVertexFistLayer.

        Vector prov(rhs.nlin());
        for (size_t s=0;s<n_dipoles;++s) {
            PROGRESSBAR(s, n_dipoles);
            const Vect3 r(dipoles(s,0),dipoles(s,1),dipoles(s,2));
            const Vect3 q(dipoles(s,3),dipoles(s,4),dipoles(s,5));

            const unsigned domainID = (dipoles_in_cortex) ? 0 : geo.getDomain(r);
            const double   sigma    = geo.sigma_in(domainID);

            unsigned istart = 0;
            prov.set(0);
            if (domainID!=0) {
                for (unsigned id=0;id<(domainID-1);++id)
                    istart += geo.getM(id).size();

                //  Treat the internal surface.

                const int nVertexLayer = geo.getM(domainID-1).nbPts();
                const int nFacesLayer = geo.getM(domainID-1).nbTrgs();
                operatorDipolePotDer(r,q,geo.getM(domainID-1),prov,istart,gauss_order,adapt_rhs);

                for(unsigned i=istart;i<istart+nVertexLayer;++i)
                    prov(i) *= -K;

                operatorDipolePot(r,q,geo.getM(domainID-1),prov,istart+nVertexLayer,gauss_order,adapt_rhs);

                for(unsigned i = istart+nVertexLayer;i<istart+nVertexLayer+nFacesLayer;++i)
                    prov(i) *= (K/sigma);

                istart += nVertexLayer+nFacesLayer;
            }

            //  Treat the external surface.

            const int nVertexLayer = geo.getM(domainID).nbPts();
            const int nFacesLayer = geo.getM(domainID).nbTrgs();

            //  Block is nVertexLayer.

            operatorDipolePotDer(r, q, geo.getM(domainID), prov, istart, gauss_order, adapt_rhs);

            for (unsigned i=istart;i<istart+nVertexLayer;++i)
                prov(i) *= K;

            //  Block is nFaceLayer.

            if (geo.nb()>(domainID+1)) {
                operatorDipolePot(r,q,geo.getM(domainID),prov,istart+nVertexLayer,gauss_order,adapt_rhs);
                for(unsigned i=istart+nVertexLayer;i<istart+nVertexLayer+nFacesLayer;++i)
                    prov(i) *= (-K/sigma);
            }
            rhs.setcol(s,prov);
        }
    }

    DipSourceMat::DipSourceMat(const Geometry& geo,const Matrix& dipoles,const int gauss_order,
                                const bool adapt_rhs,const bool dipoles_in_cortex)
    {
        assemble_DipSourceMat(*this,geo,dipoles,gauss_order,adapt_rhs,dipoles_in_cortex);
    }

    void assemble_EITSourceMat(Matrix& mat,const Geometry& geo,const Matrix& positions,const int gauss_order) {

        //  A Matrix to be applied to the scalp-injected current to obtain the Source Term of the EIT foward problem.

        const int newsize = geo.size()-geo.getM(geo.nb()-1).nbTrgs();
        mat = Matrix(newsize,positions.nlin());

        //  transmat = a big SymMatrix of which mat = part of its transpose.

        SymMatrix transmat(geo.size());
        transmat.set(0.0);

        int offset0,offset1,offset2,offset3,offset4;

        const double K = 1.0/(4*M_PI);

        for(int c=0,offset=0;c<geo.nb()-1;++c) {
            offset0 = offset;
            offset1 = offset0+geo.getM(c).nbPts();
            offset2 = offset1+geo.getM(c).nbTrgs();
            offset3 = offset2+geo.getM(c+1).nbPts();
            offset4 = offset3+geo.getM(c+1).nbTrgs();
            offset  = offset2;
        }
        const int c = geo.nb()-2;

        //  Compute S.

        operatorS(geo.getM(c+1),geo.getM(c),transmat,offset3,offset1,gauss_order);
        mult(transmat,offset3,offset1,offset4,offset2,K/geo.sigma_in(c+1));

        //  First compute D, then it will be transposed.

        operatorD(geo.getM(c+1),geo.getM(c),transmat,offset3,offset0,gauss_order);
        mult(transmat,offset3,offset0,offset4,offset1,-K);
        operatorD(geo.getM(c+1),geo.getM(c+1),transmat,offset3,offset2,gauss_order);
        mult(transmat,offset3,offset2,offset4,offset3,-2.0*K);
        operatorP1P0(geo.getM(c+1),transmat,offset3,offset2);
        mult(transmat,offset3,offset2,offset4,offset3,-0.5);

        //  Extracting the transpose of the last block of lines of transmat.
        //  Transposing the Matrix.

        for (unsigned ielec=0;ielec<positions.nlin();++ielec) {
            const Vect3 current_position(positions(ielec,0),positions(ielec,1),positions(ielec,2));
            Vect3 current_alphas; //not used here
            int current_nearest_triangle; //    To hold the index of the closest triangle to electrode.
            dist_point_mesh(current_position,geo.getM(geo.nb()-1),current_alphas,current_nearest_triangle);
            const double inv_area = 1.0/geo.getM(geo.nb()-1).triangle(current_nearest_triangle).getArea();
            for (int i=0;i<newsize;++i)
                mat(i,ielec) = transmat(newsize+current_nearest_triangle,i)*inv_area;
        }
    }

    EITSourceMat::EITSourceMat(const Geometry& geo,Sensors& electrodes,const int gauss_order) {
        assemble_EITSourceMat(*this, geo, electrodes.getPositions(), gauss_order);
    }

    void assemble_DipSource2InternalPotMat(Matrix& mat,const Geometry& geo,const Matrix& dipoles,
                                           const Matrix& points,const bool dipoles_in_cortex)
    {
        // Points with one more column for the index of the domain they belong
        Matrix pointsLabelled(points.nlin(),4);
        for (unsigned i=0;i<points.nlin();++i) {
            pointsLabelled(i,3) = geo.getDomain(Vect3(points(i,0),points(i,1),points(i,2)));
            for (int j=0;j<3;++j)
                pointsLabelled(i, j) = points(i,j);
        }
        const double K = 1.0/(4.0*M_PI);
        mat = Matrix(points.nlin(),dipoles.nlin());
        mat.set(0.0);

        for (unsigned iDIP=0;iDIP<dipoles.nlin();++iDIP) {
            const Vect3 r0(dipoles(iDIP,0),dipoles(iDIP,1),dipoles(iDIP,2));
            const Vect3  q(dipoles(iDIP,3),dipoles(iDIP,4),dipoles(iDIP,5));

            const unsigned domainID = (dipoles_in_cortex) ? 0 : geo.getDomain(r0);
            const double   sigma    = geo.sigma_in(domainID);

            static analyticDipPot anaDP;
            anaDP.init(q,r0);
            for (unsigned iPTS=0;iPTS<points.nlin();++iPTS) {
                if ((pointsLabelled(iPTS,3))==domainID) {
                    const Vect3 r(points(iPTS,0),points(iPTS,1),points(iPTS,2));
                    mat(iPTS, iDIP) = K/sigma*anaDP.f(r);
                }
            }
        }
    }

    DipSource2InternalPotMat::DipSource2InternalPotMat(const Geometry& geo,const Matrix& dipoles,
                                                       const Matrix& points,const bool dipoles_in_cortex)
    {
        assemble_DipSource2InternalPotMat(*this,geo,dipoles,points,dipoles_in_cortex);
    }

}
