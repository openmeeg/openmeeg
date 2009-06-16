/* FILE: $Id$ */

/*
Project Name : OpenMEEG

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

#include "operators.h"
#include "integrator.h"
#include "analytics.h"

namespace OpenMEEG {

    void operatorDinternal(const Mesh &m,Matrix &mat,const int offsetJ,const Matrix &points)
    {
        std::cout<<"INTERNAL OPERATOR D..."<<std::endl;
        for(size_t i=0;i<points.nlin();i++)  {
            Vect3 pt(points(i,0),points(i,1),points(i,2));
            for(int j=offsetJ;j<offsetJ+m.nbTrgs();j++){
                _operatorDinternal(i,j-offsetJ,m,mat,offsetJ,pt);
            }
        }
    }

    void operatorSinternal(const Mesh &m,Matrix &mat,const int offsetJ,const Matrix &points)
    {
        std::cout<<"INTERNAL OPERATOR S..."<<std::endl;
        for(size_t i=0;i<points.nlin();i++) {
            Vect3 pt(points(i,0),points(i,1),points(i,2));
            for(int j=offsetJ;j<offsetJ+m.nbTrgs();j++)
            {
                mat(i,j)=_operatorSinternal(j-offsetJ,m,pt);
            }
        }
    }

    // general routine for applying _operatorFerguson (see this function for further comments)
    // to an entire mesh, and storing coordinates of the output in a Matrix.
    void operatorFerguson(const Vect3& x, const Mesh &m, Matrix &mat, int offsetI, int offsetJ)
    {
        #ifdef USE_OMP
        #pragma omp parallel for
        #endif
        for(int j=offsetJ;j<offsetJ+m.nbPts();j++)
        {
            Vect3 v = _operatorFerguson(x,j-offsetJ,m);
            mat(offsetI+0,j) = v.x();
            mat(offsetI+1,j) = v.y();
            mat(offsetI+2,j) = v.z();
        }
    }

    void operatorDipolePotDer(const Vect3 &r0,const Vect3 &q,const Mesh &inner_layer,Vector &rhs,const int offsetIdx,const int GaussOrder)
    {
        static analyticDipPotDer anaDPD;

    #ifdef ADAPT_RHS
        AdaptiveIntegrator<Vect3,analyticDipPotDer> gauss(0.001);
    #else
        static Integrator<Vect3,analyticDipPotDer> gauss;
    #endif //ADAPT_RHS
        gauss.setOrder(GaussOrder);
        #ifdef USE_OMP
        #pragma omp parallel for private(anaDPD)
        #endif
        for(int i=0;i<inner_layer.nbTrgs();i++)
        {
            anaDPD.init(inner_layer,i,q,r0);
            Vect3 v=gauss.integrate(anaDPD,inner_layer.getTrg(i),inner_layer);
            #ifdef USE_OMP
            #pragma omp critical
            #endif
            {
            rhs(inner_layer.getTrg(i-offsetIdx).s1()+offsetIdx)+=v(0);
            rhs(inner_layer.getTrg(i-offsetIdx).s2()+offsetIdx)+=v(1);
            rhs(inner_layer.getTrg(i-offsetIdx).s3()+offsetIdx)+=v(2);
            }
        }
    }

    void operatorDipolePot(const Vect3 &r0, const Vect3 &q, const Mesh &inner_layer, Vector &rhs,const int offsetIdx,const int GaussOrder)
    {
        static analyticDipPot anaDP;

        anaDP.init(q,r0);
    #ifdef ADAPT_RHS
        AdaptiveIntegrator<double,analyticDipPot> gauss(0.001);
    #else
        static Integrator<double,analyticDipPot> gauss;
    #endif
        gauss.setOrder(GaussOrder);
        #ifdef USE_OMP
        #pragma omp parallel for
        #endif
        for(int i=offsetIdx;i<offsetIdx+inner_layer.nbTrgs();i++)
        {
            double d = gauss.integrate(anaDP,inner_layer.getTrg(i-offsetIdx),inner_layer);
            #ifdef USE_OMP
            #pragma omp critical
            #endif
            rhs(i) += d;
        }
    }
}

