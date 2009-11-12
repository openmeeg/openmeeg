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

/*! \file
    \brief file containing the integral operators
*/
#ifndef OPENMEEG_OPERATORS_H
#define OPENMEEG_OPERATORS_H

#include <iostream>

#include "vector.h"
#include "matrix.h"
#include "symmatrix.h"
#include "geometry.h"
#include "integrator.h"
#include "analytics.h"

namespace OpenMEEG {

    #define OPTIMIZED_OPERATOR_N

    #ifdef USE_OMP
    #undef OPTIMIZED_OPERATOR_D
    #else
    #define OPTIMIZED_OPERATOR_D
    #endif

    //#define ADAPT_LHS

    // T can be a Matrix or SymMatrix
    template<class T>
    void operatorN(const Mesh &m1,const Mesh &m2,T &mat, const int offsetI,const int offsetJ,const int,const int IopS=0,const int JopS=0);
    template<class T>
    void operatorS(const Mesh &m1,const Mesh &m2,T &mat, const int offsetI,const int offsetJ,const int);
    template<class T>
    void operatorD(const Mesh &m1,const Mesh &m2,T &mat, const int offsetI,const int offsetJ,const int);
    template<class T>
    void operatorP1P0(const Mesh &,T &mat,const int offsetI,const int offsetJ);

    void operatorSinternal(const Mesh &, Matrix &,const int,const Matrix &);
    void operatorDinternal(const Mesh &, Matrix &,const int,const Matrix &);
    void operatorFerguson(const Vect3& x, const Mesh &m1, Matrix &mat, int offsetI, int offsetJ);
    void operatorDipolePotDer(const Vect3 &r0, const Vect3 &q,const Mesh &inner_layer, Vector &rhs, int offsetIdx,const int,const bool);
    void operatorDipolePot(const Vect3 &r0, const Vect3 &q,const Mesh &inner_layer, Vector &rhs, int offsetIdx,const int,const bool);

    inline void mult( SymMatrix &mat, int Istart, int Jstart, int Istop, int Jstop, double coeff)
    {
        //If the upper left corner of the block is on the diagonal line of the Matrix
        //Only the half of the block has to be treated because of the symmetric storage
        if(Istart!=Jstart)
            for(int i=Istart;i<Istop;i++)
                #ifdef USE_OMP
                #pragma omp parallel for
                #endif
                for(int j=Jstart;j<Jstop;j++)
                    mat(i,j)*=coeff;
        else
            for(int i=Istart;i<Istop;i++)
                #ifdef USE_OMP
                #pragma omp parallel for
                #endif
                for(int j=Jstart;j<=i;j++)
                    mat(i,j)*=coeff;
    }

    inline void mult( Matrix &mat, int Istart, int Jstart, int Istop, int Jstop, double coeff)
    {
        //If the upper left corner of the block is on the diagonal line of the Matrix
        //Only the half of the block has to be treated because of the symmetric storage
        for(int i=Istart;i<Istop;i++)
            #ifdef USE_OMP
            #pragma omp parallel for
            #endif
            for(int j=Jstart;j<Jstop;j++)
                mat(i,j)*=coeff;
    }

    #ifndef OPTIMIZED_OPERATOR_D
    inline double _operatorD(const int nT1,const int nP2,const Mesh &m1,const Mesh &m2,const int GaussOrder)    
    {
        // consider varying order of quadrature with the distance between T1 and T2
        const Triangle &T1=m1.getTrg(nT1);
        #ifdef USE_OMP
            analyticD analyD;
        #else
            static analyticD analyD;
        #endif
    #ifdef ADAPT_LHS
        #ifdef USE_OMP
            AdaptiveIntegrator<double,analyticD> gauss(0.005);
        #else
            static AdaptiveIntegrator<double,analyticD> gauss(0.005);
        #endif
        gauss.setOrder(GaussOrder);
    #else
        #ifdef USE_OMP
            Integrator<double,analyticD> gauss(GaussOrder);
        #else
            static Integrator<double,analyticD> gauss(GaussOrder);
        #endif
    #endif //ADAPT_LHS

        double total = 0;

        const intSet& Tadj = m2.getTrianglesForPoint(nP2); // loop on triangles of which nP2 is a vertex
        for(intSet::const_iterator it = Tadj.begin(); it != Tadj.end(); ++it)
        {
            analyD.init( m2, *it, nP2);
            total += gauss.integrate(analyD,T1,m1);
        }

        return total;
    }
    #else

    template<class T>
    inline void _operatorD(const int nT1,const int nT2,const Mesh &m1,const Mesh &m2,T &mat,const int offsetI,const int offsetJ,const int GaussOrder)
    {
        //this version of _operatorD add in the Matrix the contribution of T2 on T1
        // for all the P1 functions it gets involved

        // consider varying order of quadrature with the distance between T1 and T2
        const Triangle &T1 = m1.getTrg(nT1);
        const Triangle &T2 = m2.getTrg(nT2);

        #ifdef USE_OMP
            analyticD3 analyD;
        #else
            static analyticD3 analyD;
        #endif

        analyD.init( m2, nT2);
    #ifdef ADAPT_LHS
        AdaptiveIntegrator<Vect3,analyticD3> gauss(0.005);
        gauss.setOrder(GaussOrder);
        Vect3 total = gauss.integrate(analyD,T1,m1);
    #else
        #ifdef USE_OMP
            Integrator<Vect3,analyticD3> gauss(GaussOrder);
        #else
            static Integrator<Vect3,analyticD3> gauss(GaussOrder);
        #endif

        Vect3 total = gauss.integrate(analyD,T1,m1);
    #endif //ADAPT_LHS

        mat(offsetI+nT1,offsetJ+((Triangle)T2)[0]) += total.x();
        mat(offsetI+nT1,offsetJ+((Triangle)T2)[1]) += total.y();
        mat(offsetI+nT1,offsetJ+((Triangle)T2)[2]) += total.z();
    }
    #endif //OPTIMIZED_OPERATOR_D


    inline void _operatorDinternal(const int nT,const int nT2, const Mesh &m, Matrix  &mat,  const int offsetJ,const Vect3 pt)
    {
        const Triangle &T2=m.getTrg(nT2); 
        static analyticD3 analyD;
      
        //vect3 total(0,0,0);
      
        analyD.init( m, nT2);

            Vect3 total=analyD.f(pt);
            
        mat(nT,offsetJ+((Triangle)T2).som(1))+=total.x();
        mat(nT,offsetJ+((Triangle)T2).som(2))+=total.y();
        mat(nT,offsetJ+((Triangle)T2).som(3))+=total.z();
     }

    inline double _operatorS(const int nT1,const int nT2,const Mesh &m1,const Mesh &m2,const int GaussOrder)
    {
        const Triangle &T1=m1.getTrg(nT1);
        const Triangle &T2=m2.getTrg(nT2);

        #ifdef USE_OMP
            Triangle *oldT=0;
            analyticS analyS;
        #else
            static Triangle *oldT=0;
            static analyticS analyS;
        #endif

        if(oldT != &T1) // a few computations are needed only when changing triangle T1
        {
            oldT=(Triangle*)&T1;
            analyS.init(nT1,m1);
        }
    #ifdef ADAPT_LHS
        AdaptiveIntegrator<double,analyticS> gauss(0.005);
        gauss.setOrder(GaussOrder);
        return gauss.integrate(analyS,T2,m2);
    #else
        #ifdef USE_OMP
            Integrator<double,analyticS> gauss;
        #else
            static Integrator<double,analyticS> gauss;
        #endif
        gauss.setOrder(GaussOrder);
        return gauss.integrate(analyS,T2,m2);
    #endif //ADAPT_LHS
    }

    inline double _operatorSinternal(const int nT, const Mesh &m,const Vect3 pt)
    {
    static analyticS analyS;
    analyS.init(nT,m);
    return analyS.f(pt);
    }

    template<class T> 
    inline double _operatorN(const int nP1,const int nP2,const Mesh &m1,const Mesh &m2,const int GaussOrder,const int IopS,const int JopS,const T &mat)
    {
        const Vect3 P1=m1.getPt(nP1);
        const Vect3 P2=m2.getPt(nP2);

        double Iqr,Aqr;
        double result=0.0;

        const intSet& trgs1 = m1.getTrianglesForPoint(nP1);
        const intSet& trgs2 = m2.getTrianglesForPoint(nP2);
        for(intSet::const_iterator it1 = trgs1.begin(); it1 != trgs1.end(); ++it1)
            for(intSet::const_iterator it2 = trgs2.begin(); it2 != trgs2.end(); ++it2)
            {
                const Triangle& T1 = m1.getTrg(*it1);
                const Triangle& T2 = m2.getTrg(*it2);

                // A1 , B1 , A2, B2 are the two opposite vertices to P1 and P2 (triangles A1,B1,P1 and A2,B2,P2)
                if(IopS!=0 || JopS!=0) Iqr=mat(IopS + *it1,JopS + *it2); else Iqr=_operatorS(*it1,*it2,m1,m2,GaussOrder);
                int nP1T=T1.contains(nP1);    //index of P1 in current triangle of mesh m1
                int nP2T=T2.contains(nP2);    //index of P2 in current triangle of mesh m2
    #ifndef OPTIMIZED_OPERATOR_N
                Vect3 A1=m1.getPt(((Triangle)T1).next(nP1T));
                Vect3 B1=m1.getPt(((Triangle)T1).prev(nP1T));
                Vect3 A2=m2.getPt(((Triangle)T2).next(nP2T));
                Vect3 B2=m2.getPt(((Triangle)T2).prev(nP2T));
                Vect3 A1B1=B1-A1;
                Vect3 A2B2=B2-A2;
                Vect3 A1P1=P1-A1;
                Vect3 A2P2=P2-A2;
                double coef1=A1P1*A1B1*1.0/A1B1.norm2();
                double coef2=A2P2*A2B2*1.0/A2B2.norm2();
                Vect3 aq=P1-(A1+A1B1*coef1);
                Vect3 br=P2-(A2+A2B2*coef2);
                aq*=(1.0/aq.norm2());
                br*=(1.0/br.norm2());

                Aqr=-0.25/T1.getArea()/T2.getArea()*( (aq^T1.normal())*(br^T2.normal()) );
    #else
                Vect3 CB1=m1.getPt(((Triangle)T1).next(nP1T))-m1.getPt(((Triangle)T1).prev(nP1T));
                Vect3 CB2=m2.getPt(((Triangle)T2).next(nP2T))-m2.getPt(((Triangle)T2).prev(nP2T));

                Aqr=-0.25/T1.getArea()/T2.getArea()*( CB1*CB2);
    #endif
                result+=Aqr*Iqr;
            }

            return result;
    }

    inline double _operatorP1P0( int nT2, int nP1, const Mesh &m)
    {
        const Triangle &T2=m.getTrg(nT2);
        if(T2.contains(nP1)== 0) {
            return 0;
        }
        else return T2.getArea()/3.;
    }

    template<class T>
    void operatorN(const Mesh &m1,const Mesh &m2,T &mat,const int offsetI,const int offsetJ,const int GaussOrder,const int IopS,const int JopS)
    {
        // This function has the following arguments:
        //    One geometry
        //    the indices of the treated layers I and J
        //    the storage Matrix for the result
        //    the upper left corner of the submatrix to be written is the Matrix
        //  the upper left corner of the corresponding S block

        std::cout<<"OPERATOR N... (arg : mesh m1, mesh m2)"<<std::endl;

        if(&m1==&m2) {
            for(int i=offsetI;i<offsetI+m1.nbPts();i++) {
                progressbar(i-offsetI,m1.nbPts());
                #ifdef USE_OMP
                #pragma omp parallel for
                #endif
                for(int j=i;j<offsetJ+m2.nbPts();j++)
                {
                    mat(i,j)=_operatorN(i-offsetI,j-offsetJ,m1,m2,GaussOrder,IopS,JopS,mat);
                }
            }
        } else {
            for(int i=offsetI;i<offsetI+m1.nbPts();i++){
                progressbar(i-offsetI,m1.nbPts());
                #ifdef USE_OMP
                #pragma omp parallel for
                #endif
                for(int j=offsetJ;j<offsetJ+m2.nbPts();j++)
                {
                    mat(i,j)=_operatorN(i-offsetI,j-offsetJ,m1,m2,GaussOrder,IopS,JopS,mat);
                }
            }
        }
    }

    template<class T> 
    void operatorS(const Mesh &m1,const Mesh &m2,T &mat,const int offsetI,const int offsetJ,const int GaussOrder)
    {
        // This function has the following arguments:
        //    One geometry
        //    the indices of the treated layers I and J
        //    the storage Matrix for the result
        //    the upper left corner of the submatrix to be written is the Matrix

        std::cout<<"OPERATOR S... (arg : mesh m1, mesh m2)"<<std::endl;

        // The operator S is given by Sij=\Int G*PSI(I,i)*Psi(J,j) with PSI(A,a) is a P0 test function on layer A and triangle a
        if(&m1==&m2)
            for(int i=offsetI;i<offsetI+m1.nbTrgs();i++) {
                progressbar(i-offsetI,m1.nbTrgs());
                #ifdef USE_OMP
                #pragma omp parallel for
                #endif
                for(int j=i;j<offsetJ+m2.nbTrgs();j++)
                {
                    mat(i,j)=_operatorS(i-offsetI,j-offsetJ,m1,m2,GaussOrder);
                }
            }
        else
        {
            for(int i=offsetI;i<offsetI+m1.nbTrgs();i++) {
                progressbar(i-offsetI,m1.nbTrgs());
                #ifdef USE_OMP
                #pragma omp parallel for
                #endif
                for(int j=offsetJ;j<offsetJ+m2.nbTrgs();j++)
                {
                    mat(i,j)=_operatorS(i-offsetI,j-offsetJ,m1,m2,GaussOrder);
                }
            }
        }
    }

    #ifndef OPTIMIZED_OPERATOR_D

    template<class T>
    void operatorD(const Mesh &m1,const Mesh &m2,T &mat,const int offsetI,const int offsetJ,const int GaussOrder)
    // This function (NON OPTIMIZED VERSION) has the following arguments:
    //    One geometry
    //    the indices of the treated layers I and J
    //    the storage Matrix for the result
    //    the upper left corner of the submatrix to be written is the Matrix
    {
        std::cout<<"OPERATOR D... (arg : mesh m1, mesh m2)"<<std::endl;

        for(int i=offsetI;i<offsetI+m1.nbTrgs();i++) {
            progressbar(i-offsetI,m1.nbTrgs());
            #ifdef USE_OMP
            #pragma omp parallel for
            #endif
            for(int j=offsetJ;j<offsetJ+m2.nbPts();j++)
            {
                // P1 functions are tested thus looping on vertices
                mat(i,j)=_operatorD(i-offsetI,j-offsetJ,m1,m2,GaussOrder);
            }
        }
    }

    #else // OPTIMIZED_OPERATOR_D

    template<class T>
    void operatorD(const Mesh &m1,const Mesh &m2,T &mat,const int offsetI,const int offsetJ,const int GaussOrder)
    {
        // This function (OPTIMIZED VERSION) has the following arguments:
        //    One geometry
        //    the indices of the treated layers I and J
        //    the storage Matrix for the result
        //    the upper left corner of the submatrix to be written is the Matrix

        std::cout<<"OPERATOR D (Optimized) ... (arg : mesh m1, mesh m2)"<<std::endl;

        for(int i=offsetI;i<offsetI+m1.nbTrgs();i++) {
            progressbar(i-offsetI,m1.nbTrgs());
            for(int j=offsetJ;j<offsetJ+m2.nbTrgs();j++)
            {
                //In this version of the funtcion, in order to skip multiple computations of the same quantities
                //    loops are run over the triangles but the Matrix cannot be filled in this function anymore
                //    That's why the filling is done is function _operatorD
                _operatorD(i-offsetI,j-offsetJ,m1,m2,mat,offsetI,offsetJ,GaussOrder);
            }
        }
    }

    #endif // OPTIMIZED_OPERATOR_D

    template<class T>
    void operatorP1P0(const Mesh &m,T &mat,const int offsetI,const int offsetJ)
    {
        // This time mat(i,j)+= ... the Matrix is incremented by the P1P0 operator
        std::cout<<"OPERATOR P1P0..."<<std::endl;
        for(int i=offsetI;i<offsetI+m.nbTrgs();i++)
            for(int j=offsetJ;j<offsetJ+m.nbPts();j++)
            {
                mat(i,j)+=_operatorP1P0(i-offsetI,j-offsetJ,m);
            }
    }


    inline Vect3 _operatorFerguson(const Vect3 x,const int nP1,const Mesh &m1)
    {
        const Vect3 P1=m1.getPt(nP1);

        double opS;
        Vect3  v;

        #ifdef USE_OMP
            Vect3 result;
            analyticS analyS;
        #else
            static Vect3 result;
            static analyticS analyS;
        #endif

        result.x()=0;
        result.y()=0;
        result.z()=0;

        //loop over triangles of which P1 is a vertex
        const intSet& trgs1 = m1.getTrianglesForPoint(nP1);
        for(intSet::const_iterator it = trgs1.begin(); it != trgs1.end(); ++it)
        {
            const Triangle& T1=m1.getTrg(*it);

            // A1 , B1  are the two opposite vertices to P1 (triangles A1,B1,P1)
            int nP1T = T1.contains(nP1);    //index of P1 in current triangle of mesh m1

            Vect3 A1 = m1.getPt(T1.next(nP1T));
            Vect3 B1 = m1.getPt(T1.prev(nP1T));
            Vect3 A1B1 = B1-A1;    // actually, B1A1 is needed
            v = A1B1*(-0.5/T1.getArea());

            analyS.init(P1,A1,B1);
            opS = analyS.f(x);

            result += (v*opS);
        }

        return result;
    }
}

#endif  //! OPENMEEG_OPERATORS_H
