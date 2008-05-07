/* FILE: $Id$ */

/*
Project Name : OpenMEEG

author            : $Author$
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

#ifndef H_integrator
#define H_integrator

#include "vect3.h"
#include "triangle.h"
#include "mesh3.h"
#include <cmath>
#include <iostream>

// light class containing d Vect3
template <int d>
class Vect3array
{
private:
    Vect3 t[d];

public:
    Vect3array() {};
    inline Vect3array(double x) {
        for (int i=0;i<d;i++)
            t[i]=Vect3(x);
    }
    inline Vect3array<d> operator*(double x) const {
        Vect3array<d> r;
        for (int i=0;i<d;i++)
            r.t[i]=t[i]*x;
        return r;
    }
    inline Vect3 operator()(int i) const { return t[i]; }
    inline Vect3& operator()(int i) { return t[i]; }
};

template <int d>
inline void multadd (Vect3array<d> &target, const double scale,  const Vect3array<d> &incr)
{
    for (int i=0;i<d;i++) {
        target(i) = target(i) + scale*incr(i);
    }
}

inline void multadd (double &target, const double scale, const double incr)
{
    target += scale*incr;
}

inline void multadd (Vect3 &target, const double scale,  const Vect3 &incr)
{
    target = target + scale*incr;
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

template<class T,class I>
class integrator
{
private:
    int order;

public:
    inline integrator() {setOrder(3);}
    inline integrator(int ord) {setOrder(ord);}
    inline ~integrator() {}
    inline void setOrder(int n)
    {
        if(n>=0 && n<4) order=n;
        else {std::cout<<"Unavalaible Gauss Order: "<<n<<std::endl; order = (n<1)?order=1:order;}
    }

    inline T integrate ( const I &fc, const Triangle& Trg ,const Mesh& M)
    {
        Vect3 sommets[3]={M.getPt(Trg.s1()),M.getPt(Trg.s2()),M.getPt(Trg.s3())};
        return triangle_integration(fc,sommets);
    }
    inline T triangle_integration( const I &fc, Vect3 *vertices)
    {// compute double area of triangle defined by vertices
        Vect3 crossprod=(vertices[1]-vertices[0])^(vertices[2]-vertices[0]);
        double S = crossprod.norme();
        T result = 0;
        static Vect3 zero(0.0,0.0,0.0);
        int i;
        for(i=0;i<nbPts[order];i++)
        {
            Vect3 v=zero;
            int j;
            for(j=0;j<3;j++) {
                v.multadd(cordBars[order][i][j],vertices[j]);
            }
            multadd(result,cordBars[order][i][3],fc.f(v));
        }
        return result*S;
    }
};

template<class T,class I>
class adaptive_integrator : public integrator<T,I>
{
private:
    double tolerance;
public:
    inline adaptive_integrator() : tolerance(0.0001) {}
    inline adaptive_integrator(double tol) : tolerance(tol) {}
    inline ~adaptive_integrator() {}
    inline double norme(double a) {
        return fabs(a);
    }
    inline double norme(Vect3 a) {
        return a.norme();
    }
    inline T integrate(const I &fc, const Triangle& Trg ,const Mesh& M)
    {
        int n=0;
        Vect3 vertices[3]={M.getPt(Trg.s1()),M.getPt(Trg.s2()),M.getPt(Trg.s3())};
        T I0=triangle_integration(fc,vertices);
        return adaptive_integration(fc,vertices,I0,tolerance,n);
    }
    inline T adaptive_integration(const I &fc,const Vect3 *vertices,T I0,const double tolerance,int n)
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
        T sum=I1+I2+I3+I4;
        if (norme(I0-sum)>tolerance*norme(I0)){
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

