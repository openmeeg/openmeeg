/*! \file
    \brief file containing the integral operators
*/
#ifndef H_operateurs
#define H_operateurs

#include "vecteur.h"
#include "matrice.h"
#include "symmatrice.h"
#include "geometry.h"
#include "integrateur.h"
#include "analytiques.h"

#define OPTIMIZED_OPERATOR_N
#define OPTIMIZED_OPERATOR_D

void operateurN(const Geometry &geo,const int I,const int J,const int,symmatrice &mat,const int offsetI,const int offsetJ,const int IopS,const int JopS);
void operateurS(const Geometry &geo,const int I,const int J,const int,symmatrice &mat,const int offsetI,const int offsetJ);
void operateurD(const Geometry &geo,const int I,const int J,const int,symmatrice &mat,const int offsetI,const int offsetJ);

void operateurN(const Mesh &m1,const Mesh &m2,genericMatrix &mat,const int offsetI,const int offsetJ,const int,const int IopS=0,const int JopS=0);
void operateurS(const Mesh &m1,const Mesh &m2,genericMatrix &mat,const int offsetI,const int offsetJ,const int);
void operateurD(const Mesh &m1,const Mesh &m2,genericMatrix &mat,const int offsetI,const int offsetJ,const int);

void operateurFerguson(const Vect3 x, const Mesh &m1, matrice &mat, int offsetI, int offsetJ);
void trBlock(matrice &mat,int Iinit,int Jinit,int nbI,int nbJ);
void operateurDipolePotDer(const Vect3 &r0, const Vect3 &q,const Mesh &inner_layer, vecteur &rhs, int offsetIdx,const int);
void operateurDipolePot(const Vect3 &r0, const Vect3 &q,const Mesh &inner_layer, vecteur &rhs, int offsetIdx,const int);
void operateurDipolePotDerGrad(const Vect3 &r0, const Vect3 &q,const Mesh &inner_layer, vecteur rhs[3], int offsetIdx,const int);
void operateurDipolePotGrad(const Vect3 &r0, const Vect3 &q,const Mesh &inner_layer, vecteur rhs[3], int offsetIdx,const int);

inline void mult( symmatrice &mat, int Istart, int Jstart, int Istop, int Jstop, double coeff)
{
    //If the upper left corner of the block is on the diagonal line of the matrix
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
            for(int j=Jstart;j<=Jstart+(i-Istart);j++)
                mat(i,j)*=coeff;

}

inline void mult2( matrice &mat, int Istart, int Jstart, int Istop, int Jstop, double coeff)
{
    //If the upper left corner of the block is on the diagonal line of the matrix
    //Only the half of the block has to be treated because of the symmetric storage
    for(int i=Istart;i<Istop;i++)
        #ifdef USE_OMP
        #pragma omp parallel for
        #endif
        for(int j=Jstart;j<Jstop;j++)
            mat(i,j)*=coeff;
}

#ifndef OPTIMIZED_OPERATOR_D
inline double _operateurD(const int nT1,const int nP2,const int GaussOrder,const Mesh &m1,const Mesh &m2)
{
    // consider varying order of quadrature with the distance between T1 and T2
    const Triangle &T1=m1.getTrg(nT1);

    static integrateur<double> gauss(GaussOrder);
    static analytiqueD analyD;

    double total = 0;
    static int Tadj[128];            // triangles of which P2 is a vertex
    int nTadj = m2.elem(nP2,Tadj);

    for(int k=0;k<nTadj;k++)
    {
        analyD.init( m2, Tadj[k], nP2);
        total += gauss.integre(analyD,T1,m1);
    }

    return total;
}
#else
inline void _operateurD(const int nT1,const int nT2,const int GaussOrder,const Mesh &m1,const Mesh &m2,genericMatrix &mat,const int offsetI,const int offsetJ)
{
    //this version of _operateurD add in the matrix the contribution of T2 on T1
    // for all the P1 functions it gets involved

    // consider varying order of quadrature with the distance between T1 and T2
    const Triangle &T1=m1.getTrg(nT1);
    const Triangle &T2=m2.getTrg(nT2);

    static integrateur<Vect3> gauss(GaussOrder);
    static analytiqueD3 analyD;

    analyD.init( m2, nT2);
    Vect3 total = gauss.integre(analyD,T1,m1);
    mat(offsetI+nT1,offsetJ+((Triangle)T2)[0]) += total.x();
    mat(offsetI+nT1,offsetJ+((Triangle)T2)[1]) += total.y();
    mat(offsetI+nT1,offsetJ+((Triangle)T2)[2]) += total.z();
}
#endif

inline double _operateurS(const int nT1,const int nT2,const int GaussOrder,const Mesh &m1,const Mesh &m2)
{
    const Triangle &T1=m1.getTrg(nT1);
    const Triangle &T2=m2.getTrg(nT2);

    static Triangle *oldT=0;
    static integrateur<double> gauss;
    static analytiqueS analyS;
    gauss.setOrdre(GaussOrder);

    if(oldT != &T1)    // a few computations are needed only when changing triangle T1
    {
        oldT=(Triangle*)&T1;
        analyS.init(nT1,m1);
    }

    return gauss.integre(analyS,T2,m2);
}

inline double _operateurN(const int nP1,const int nP2,const int GaussOrder,const Mesh &m1,const Mesh &m2,const int IopS,const int JopS,genericMatrix &mat)
{
    int trgs1[128],trgs2[128];
    const Vect3 P1=m1.getPt(nP1);
    const Vect3 P2=m2.getPt(nP2);

    int nik=m1.elem(nP1,trgs1);    //number of triangles of which P1 is a vertex
    int njl=m2.elem(nP2,trgs2);    //number of triangles of which P2 is a vertex

    double Iqr,Aqr;
    double result=0.0;

    for(int q=0;q<nik;q++)        //loop over triangles of which P1 is a vertex
        for(int r=0;r<njl;r++)    //loop over triangles of which P2 is a vertex
        {
            const Triangle& T1=m1.getTrg(trgs1[q]);
            const Triangle& T2=m2.getTrg(trgs2[r]);

            // A1 , B1 , A2, B2 are the two opposite vertices to P1 and P2 (triangles A1,B1,P1 and A2,B2,P2)
            if(IopS!=0 || JopS!=0) Iqr=mat(IopS+trgs1[q],JopS+trgs2[r]); else Iqr=_operateurS(trgs1[q],trgs2[r],GaussOrder,m1,m2);
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
            double coef1=A1P1*A1B1*1.0/A1B1.norme2();
            double coef2=A2P2*A2B2*1.0/A2B2.norme2();
            Vect3 aq=P1-(A1+A1B1*coef1);
            Vect3 br=P2-(A2+A2B2*coef2);
            aq*=(1.0/aq.norme2());
            br*=(1.0/br.norme2());

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


//calcultates the S at point x integrated over all the triagles having nP1 as a vertice.
inline Vect3& _operateurFerguson(const Vect3 x,const int nP1,const Mesh &m1)
{
    int trgs1[128];
    const Vect3 P1=m1.getPt(nP1);

    int nik=m1.elem(nP1,trgs1);    //number of triangles of which P1 is a vertex

    double opS;
    Vect3  v;
    static Vect3 result;
    static analytiqueS analyS;

    result.x()=0;
    result.y()=0;
    result.z()=0;

    for(int q=0;q<nik;q++)        //loop over triangles of which P1 is a vertex
    {
        const Triangle& T1=m1.getTrg(trgs1[q]);

        // A1 , B1 , A2, B2 are the two opposite vertices to P1 and P2 (triangles A1,B1,P1 and A2,B2,P2)
        int nP1T=T1.contains(nP1);    //index of P1 in current triangle of mesh m1

        Vect3 A1=m1.getPt(((Triangle)T1).next(nP1T));
        Vect3 B1=m1.getPt(((Triangle)T1).prev(nP1T));
        Vect3 A1B1=B1-A1;    // actually, B1A1 is needed
        v=A1B1*(-0.5/T1.getArea());

        analyS.init(P1,A1,B1);
        opS=analyS.f(x);

        result+=(v*opS);
    }

    return result;
}

#endif
