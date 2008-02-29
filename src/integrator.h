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

#ifndef H_integrator
#define H_integrator

#include "vect3.h"
#include "triangle.h"
#include "mesh3.h"
#include "fcontainer.h"
#include <cmath>
#include <iostream>

inline void multadd (double &cible, const double multiplicateur, const double multiplie)
{
    cible+=multiplicateur*multiplie;
}

inline void multadd (Vect3 &cible, const double multiplicateur,  const Vect3 &multiplie)
{
    cible=cible+multiplicateur*multiplie;
}

// light class containing d Vect3
template <int d> class vect3array {
    Vect3 t[d];

public:
    vect3array() {};
    inline vect3array(double x) {
        for (int i=0;i<d;i++)
            t[i]=Vect3(x);
    }
    inline vect3array<d> operator*(double x) const {
        vect3array<d> r;
        for (int i=0;i<d;i++)
            r.t[i]=t[i]*x;
        return r;
    }
    inline Vect3 operator()(int i) const { return t[i]; }
    inline Vect3& operator()(int i) { return t[i]; }
};

template <int d>
inline void multadd (vect3array<d> &cible, const double multiplicateur,  const vect3array<d> &multiplie)
{
    for (int i=0;i<d;i++) {
        cible(i)=cible(i)+multiplicateur*multiplie(i);
    }
}

// Quadrature rules are from Marc Bonnet's book: Equations integrales..., Appendix B.3

static const double cordBars[4][16][4]=
{
    //parameters for N=3
    {
        {0.166666666666667,0.166666666666667,0.666666666666667,0.166666666666667},
        {0.166666666666667,0.666666666666667,0.166666666666667,0.166666666666667},
        {0.666666666666667,0.166666666666667,0.166666666666667,0.166666666666667},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0}
    }
    ,
    // parameters for N=6
    {
        {0.445948490915965,0.445948490915965,0.10810301816807,0.111690794839005},
        {0.445948490915965,0.10810301816807,0.445948490915965,0.111690794839005},
        {0.10810301816807,0.445948490915965,0.445948490915965,0.111690794839005},
        {0.091576213509771,0.091576213509771,0.81684757298045796,0.054975871827661},
        {0.091576213509771,0.81684757298045796,0.091576213509771,0.054975871827661},
        {0.816847572980458,0.091576213509771,0.091576213509771,0.054975871827661},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0}
    }
    ,
        // parameters for N=7
    {
        {0.333333333333333,0.333333333333333,0.333333333333333,0.1125},
        {0.470142064105115,0.470142064105115,0.059715871789770,0.066197076394253},
        {0.470142064105115,0.059715871789770,0.470142064105115,0.066197076394253},
        {0.059715871789770,0.470142064105115,0.470142064105115,0.066197076394253},
        {0.101286507323456,0.101286507323456,0.79742698535308798,0.062969590272414},
        {0.101286507323456,0.7974269853530880,0.101286507323456,0.062969590272414},
        {0.797426985353088,0.101286507323456,0.101286507323456,0.062969590272414},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0},
        {0.0,0.0,0.0,0.0}
    }
    ,

        // parameters for N=16
    {
        {0.333333333333333,0.333333333333333,0.3333333333333333,0.072157803838893},
        {0.081414823414554,0.459292588292722,0.459292588292722,0.047545817133642},
        {0.459292588292722,0.081414823414554,0.459292588292722,0.047545817133642},
        {0.459292588292722,0.459292588292722,0.081414823414554,0.047545817133642},
        {0.898905543365937,0.050547228317031,0.050547228317031,0.016229248811599},
        {0.050547228317031,0.898905543365937,0.050547228317031,0.016229248811599},
        {0.050547228317031,0.050547228317031,0.898905543365937,0.016229248811599},
        {0.658861384496479,0.170569307751760,0.17056930775176099,0.051608685267359},
        {0.170569307751760,0.658861384496479,0.17056930775176099,0.051608685267359},
        {0.170569307751760,0.17056930775176099,0.658861384496479,0.051608685267359},
        {0.008394777409957,0.728492392955404,0.263112829634639,0.013615157087217},
        {0.728492392955404,0.008394777409957,0.263112829634639,0.013615157087217},
        {0.728492392955404,0.263112829634639,0.008394777409957,0.013615157087217},
        {0.008394777409957,0.263112829634639,0.728492392955404,0.013615157087217},
        {0.263112829634639,0.008394777409957,0.728492392955404,0.013615157087217},
        {0.263112829634639,0.728492392955404,0.008394777409957,0.013615157087217}
    }

}; // end of gaussTriangleParams

static const int nbPts[4]={3,6,7,16};

template<class T> class integrator
{
private:
    // ordre numéro_du_noeud x_y_z_bary
    int ordre;

public:
    inline integrator() {setOrder(3);}
    inline integrator(int ord) {setOrder(ord);}
    inline ~integrator() {}
    inline void setOrder(int n)
    {
        if(n>=0 && n<4) ordre=n;
        else {std::cout<<"Unavalaible Gauss Order: "<<n<<std::endl; ordre = (n<1)?ordre=1:ordre;}
    }

    inline T integrate ( const fContainer<T> &fc, const Triangle& Trg ,const Mesh& M)
    {
        Vect3 sommets[3]={M.getPt(Trg.s1()),M.getPt(Trg.s2()),M.getPt(Trg.s3())};
        return triangle_integration(fc,sommets);
    }
    inline T triangle_integration( const fContainer<T> &fc, Vect3 *vertices)
    {// compute double area of triangle defined by vertices
        Vect3 crossprod=(vertices[1]-vertices[0])^(vertices[2]-vertices[0]);
        double S = crossprod.norme();
        T result = 0;
        static Vect3 zero(0.0,0.0,0.0);
        int i;
        for(i=0;i<nbPts[ordre];i++)
        {
            Vect3 v=zero;
            int j;
            for(j=0;j<3;j++) {
                v.multadd(cordBars[ordre][i][j],vertices[j]);
            }
            multadd(result,cordBars[ordre][i][3],fc.f(v));
        }
        return result*S;
    }
};

template<class T> class adaptive_integrator : public integrator<T>
{
private:
    double tolerance;
public:
    inline adaptive_integrator() {setTol(0.0001);}
    inline adaptive_integrator(double tol) {setTol(tol);}
    inline ~adaptive_integrator() {}
    inline void setTol(double tol)
    {
        tolerance = tol;
    }
    inline double norme(double a) {
        return fabs(a);
    }
    inline double norme(Vect3 a) {
        return a.norme();
    }
    inline T integrate ( const fContainer<T> &fc, const Triangle& Trg ,const Mesh& M)
    {
        int n=0;
        Vect3 sommets[3]={M.getPt(Trg.s1()),M.getPt(Trg.s2()),M.getPt(Trg.s3())};
        T I0=triangle_integration(fc,sommets);
        return adaptive_integration(fc,sommets,I0,tolerance,n);
    }
    inline T adaptive_integration(const fContainer<T> &fc,const Vect3 *vertices,T I0,const double tolerance,int n)
    {
        Vect3 newpoint0(0.0,0.0,0.0);
        multadd(newpoint0,0.5,vertices[0]);
        multadd(newpoint0,0.5,vertices[1]);
        Vect3 newpoint1(0.0,0.0,0.0);
        multadd(newpoint1,0.5,vertices[1]);
        multadd(newpoint1,0.5,vertices[2]);
        Vect3 newpoint2(0.0,0.0,0.0);
        multadd(newpoint2,0.5,vertices[2]);
        multadd(newpoint2,0.5,vertices[0]);
        Vect3 vertices1[3]={vertices[0],newpoint0,newpoint2};
        Vect3 vertices2[3]={vertices[1],newpoint1,newpoint0};
        Vect3 vertices3[3]={vertices[2],newpoint2,newpoint1};
        Vect3 vertices4[3]={newpoint0,newpoint1,newpoint2};
        T I1=triangle_integration(fc,vertices1);
        T I2=triangle_integration(fc,vertices2);
        T I3=triangle_integration(fc,vertices3);
        T I4=triangle_integration(fc,vertices4);
        T somme=I1+I2+I3+I4;
        if (norme(I0-somme)>tolerance*norme(I0)){
            n=n+1;
            if (n<10) {
                I1 = adaptive_integration(fc,vertices1,I1,tolerance,n);
                I2 = adaptive_integration(fc,vertices2,I2,tolerance,n);
                I3 = adaptive_integration(fc,vertices3,I3,tolerance,n);
                I4 = adaptive_integration(fc,vertices4,I4,tolerance,n);
                I0 = I1+I2+I3+I4;
            }
        }
        return I0;
    }
};

#endif

