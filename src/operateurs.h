#ifndef H_operateurs
#define H_operateurs

#include "vecteur.h"
#include "matrice.h"
#include "symmatrice.h"
#include "geometry.h"
#include "integrateur.h"
#include "analytiques.h"
extern int GaussOrder;
#define OPTIMIZED_OPERATOR_N
#define OPTIMIZED_OPERATOR_D

void operateurN(geometry &geo,int I,int J,symmatrice &mat,int offsetI,int offsetJ,int IopS, int JopS );
void operateurN(mesh &m1,mesh &m2,genericMatrix &mat,int offsetI,int offsetJ,int IopS=0, int JopS=0 );
void operateurS(geometry &geo,int I,int J,symmatrice &mat,int offsetI,int offsetJ);
void operateurS(const mesh &m1,const mesh &m2,genericMatrix &mat,int offsetI,int offsetJ);
void operateurD(geometry &geo,int I,int J,symmatrice &mat,int offsetI,int offsetJ);
void operateurD(mesh &m1,mesh &m2,genericMatrix &mat,int offsetI,int offsetJ);
void operateurFerguson( const vect3 x, const mesh &m1, matrice &mat, int offsetI, int offsetJ);
void trBlock(matrice &mat,int Iinit,int Jinit,int nbI,int nbJ);
void operateurDipolePotDer(const vect3 &r0, const vect3 &q,const mesh &inner_layer, vecteur &rhs, int offsetIdx);
void operateurDipolePot(const vect3 &r0, const vect3 &q,const mesh &inner_layer, vecteur &rhs, int offsetIdx);
void operateurDipolePotDerGrad(const vect3 &r0, const vect3 &q,const mesh &inner_layer, vecteur rhs[3], int offsetIdx);
void operateurDipolePotGrad(const vect3 &r0, const vect3 &q,const mesh &inner_layer, vecteur rhs[3], int offsetIdx);


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
inline double _operateurD( int nT1, int nP2, const mesh &m1, const mesh &m2)
{
    // consider varying order of quadrature with the distance between T1 and T2
    const triangle &T1=m1.trngl(nT1);

    static integrateur<double> gauss(GaussOrder);
    static analytiqueD analyD;

    double total=0;
    static int Tadj[128];            // triangles of which P2 is a vertex
    int nTadj=m2.elem(nP2,Tadj);    //

    for(int k=0;k<nTadj;k++)
    {
        analyD.init( m2, Tadj[k], nP2);
        total+=gauss.integre(analyD,T1,m1);
    }

    return total;
}
#else
inline void _operateurD( int nT1, int nT2, const mesh &m1, const mesh &m2, genericMatrix &mat, int offsetI, int offsetJ)
{
    //this version of _operateurD add in the matrix the contribution of T2 on T1
    // for all the P1 functions it gets involved

    // consider varying order of quadrature with the distance between T1 and T2
    const triangle &T1=m1.trngl(nT1);
    const triangle &T2=m2.trngl(nT2);

    static integrateur<vect3> gauss(GaussOrder);
    static analytiqueD3 analyD;

    analyD.init( m2, nT2);
    vect3 total=gauss.integre(analyD,T1,m1);
    mat(offsetI+nT1,offsetJ+((triangle)T2).som(1))+=total.X();
    mat(offsetI+nT1,offsetJ+((triangle)T2).som(2))+=total.Y();
    mat(offsetI+nT1,offsetJ+((triangle)T2).som(3))+=total.Z();
}
#endif

inline double _operateurS( int nT1, int nT2, const mesh &m1, const mesh &m2)
{
    const triangle &T1=m1.trngl(nT1);
    const triangle &T2=m2.trngl(nT2);

    static triangle *oldT=0;
    static integrateur<double> gauss;
    static analytiqueS analyS;
    gauss.setOrdre(GaussOrder);

    if(oldT != &T1)    // a few computations are needed only when changing triangle T1
    {
        oldT=(triangle*)&T1;
        analyS.init(nT1,m1);
    }

    return gauss.integre(analyS,T2,m2);
}

inline double _operateurN( int nP1, int nP2, const mesh &m1, const mesh &m2, int IopS, int JopS, genericMatrix &mat)
{
    int trgs1[128],trgs2[128];
    const vect3 P1=m1.vctr(nP1);
    const vect3 P2=m2.vctr(nP2);

    int nik=m1.elem(nP1,trgs1);    //number of triangles of which P1 is a vertex
    int njl=m2.elem(nP2,trgs2);    //number of triangles of which P2 is a vertex

    double Iqr,Aqr;
    double retour=0.0;

    for(int q=0;q<nik;q++)        //loop over triangles of which P1 is a vertex
        for(int r=0;r<njl;r++)    //loop over triangles of which P2 is a vertex
        {
            const triangle& T1=m1.trngl(trgs1[q]);
            const triangle& T2=m2.trngl(trgs2[r]);

            // A1 , B1 , A2, B2 are the two opposite vertices to P1 and P2 (triangles A1,B1,P1 and A2,B2,P2)
            if(IopS!=0 || JopS!=0) Iqr=mat(IopS+trgs1[q],JopS+trgs2[r]); else Iqr=_operateurS(trgs1[q],trgs2[r],m1,m2);
            int nP1T=T1.appar(nP1);    //index of P1 in current triangle of mesh m1
            int nP2T=T2.appar(nP2);    //index of P2 in current triangle of mesh m2
#ifndef OPTIMIZED_OPERATOR_N
            vect3 A1=m1.vctr(((triangle)T1).next(nP1T));
            vect3 B1=m1.vctr(((triangle)T1).prev(nP1T));
            vect3 A2=m2.vctr(((triangle)T2).next(nP2T));
            vect3 B2=m2.vctr(((triangle)T2).prev(nP2T));
            vect3 A1B1=B1-A1;
            vect3 A2B2=B2-A2;
            vect3 A1P1=P1-A1;
            vect3 A2P2=P2-A2;
            double coef1=A1P1*A1B1*1.0/A1B1.norme2();
            double coef2=A2P2*A2B2*1.0/A2B2.norme2();
            vect3 aq=P1-(A1+A1B1*coef1);
            vect3 br=P2-(A2+A2B2*coef2);
            aq*=(1.0/aq.norme2());
            br*=(1.0/br.norme2());

            Aqr=-0.25/T1.getAire()/T2.getAire()*( (aq^T1.normale())*(br^T2.normale()) );
#else
            vect3 CB1=m1.vctr(((triangle)T1).next(nP1T))-m1.vctr(((triangle)T1).prev(nP1T));
            vect3 CB2=m2.vctr(((triangle)T2).next(nP2T))-m2.vctr(((triangle)T2).prev(nP2T));

            Aqr=-0.25/T1.getAire()/T2.getAire()*( CB1*CB2);
#endif
            retour+=Aqr*Iqr;
        }

        return retour;
}


//calcultates the S at point x integrated over all the triagles having nP1 as a vertice.
inline vect3& _operateurFerguson( const vect3 x, int nP1, const mesh &m1)
{
    int trgs1[128];
    const vect3 P1=m1.vctr(nP1);

    int nik=m1.elem(nP1,trgs1);    //number of triangles of which P1 is a vertex

    double opS;
    vect3  v;
    static vect3 retour;
    static analytiqueS analyS;

    retour._x()=0;
    retour._y()=0;
    retour._z()=0;

    for(int q=0;q<nik;q++)        //loop over triangles of which P1 is a vertex
    {
        const triangle& T1=m1.trngl(trgs1[q]);

        // A1 , B1 , A2, B2 are the two opposite vertices to P1 and P2 (triangles A1,B1,P1 and A2,B2,P2)
        int nP1T=T1.appar(nP1);    //index of P1 in current triangle of mesh m1

        vect3 A1=m1.vctr(((triangle)T1).next(nP1T));
        vect3 B1=m1.vctr(((triangle)T1).prev(nP1T));
        vect3 A1B1=B1-A1;    // actually, B1A1 is needed
        v=A1B1*(-0.5/T1.getAire());

        analyS.init(P1,A1,B1);
        opS=analyS.f(x);

        retour+=(v*opS);
    }

    return retour;
}

#endif
