#ifndef H_analytiques
#define H_analytiques

#include "fcontainer.h"
#include "mesh3.h"

class analytiqueS : public fContainer<double> {
private:
    vect3 p0,p1,p2;    //vertices of the triangle
    vect3 p2p1,p1p0,p0p2;
    vect3 nu0,nu1,nu2;
    vect3 n;
    double norme2p2p1,norme2p1p0,norme2p0p2;
    double tanTHETA0m,tanTHETA0p,tanTHETA1m,tanTHETA1p,tanTHETA2m,tanTHETA2p;

public:
    analytiqueS(){}
    ~analytiqueS(){}
    void init( int nT, const mesh &m)
    {
        // all computations needed when the first triangle of integration is changed
        triangle& T=(triangle&)m.trngl(nT);

        p0=m.vctr(T.som1());
        p1=m.vctr(T.som2());
        p2=m.vctr(T.som3());

        p1p0=p1-p0; p2p1=p2-p1; p0p2=p0-p2;
        norme2p1p0=p1p0.norme(); norme2p2p1=p2p1.norme(); norme2p0p2=p0p2.norme();

        n=T.normale();
        n*=-(1/n.norme()) ;

        nu0=(n^p1p0);
        nu1=(n^p2p1);
        nu2=(n^p0p2);
        nu0*=(1/nu0.norme()) ;
        nu1*=(1/nu1.norme()) ;
        nu2*=(1/nu2.norme()) ;
    }

    void init( const vect3 v0, const vect3 v1,const vect3 v2)
    {
        // all computations needed when the first triangle of integration is changed

        p0=v0;
        p1=v1;
        p2=v2;

        p1p0=p1-p0; p2p1=p2-p1; p0p2=p0-p2;
        norme2p1p0=p1p0.norme(); norme2p2p1=p2p1.norme(); norme2p0p2=p0p2.norme();

        n=p1p0^p0p2;
        n*=-(1/n.norme()) ;

        nu0=(n^p1p0);
        nu1=(n^p2p1);
        nu2=(n^p0p2);
        nu0*=(1/nu0.norme()) ;
        nu1*=(1/nu1.norme()) ;
        nu2*=(1/nu2.norme()) ;
    }


    inline double f(const vect3& x) const
    {
        // analytical value of the internal integral of S operator at point X
        vect3 p1x=p1-x, p2x=p2-x, p0x=p0-x ;
        double norme2p0x=p0x.norme();
        double norme2p1x=p1x.norme();
        double norme2p2x=p2x.norme();
        double alpha=(x-p0)*n ;
        double g0,g1,g2;

        if ((p0x^p1p0).norme()>.00000001)
            g0=-log(norme2p1x-p1x*p1p0*(1.0/norme2p1p0) )+log(norme2p0x-p0x*p1p0*(1.0/norme2p1p0) );
        else
            g0= fabs(log(norme2p1x)-log(norme2p0x));
        if ((p1x^p2p1).norme()>.00000001)
            g1=-log(norme2p2x-p2x*p2p1*(1.0/norme2p2p1) )+log(norme2p1x-p1x*p2p1*(1.0/norme2p2p1) );
        else
            g1= fabs(log(norme2p2x)-log(norme2p1x));
        if ((p2x^p0p2).norme()>.00000001)
            g2=-log(norme2p0x-p0x*p0p2*(1.0/norme2p0p2) )+log(norme2p2x-p2x*p0p2*(1.0/norme2p0p2) );
        else
            g2=fabs(log(norme2p0x)-log(norme2p2x));
        return ((p0x*nu0)*g0+(p1x*nu1)*g1+(p2x*nu2)*g2)-alpha*x.solangl(p0,p1,p2);
    }
};

class analytiqueD : public fContainer<double> {
private:
    vect3 v1,v2,v3;
    int i;
    double aire;
public:
    analytiqueD(){}
    ~analytiqueD(){}
    inline void init( const mesh& m1, const int trg, const int noeud)
    {
        v1=m1.vctr(m1.trngl(trg).id1());
        v2=m1.vctr(m1.trngl(trg).id2());
        v3=m1.vctr(m1.trngl(trg).id3());
        i=m1.trngl(trg).appar(noeud);
        aire=m1.trngl(trg).getAire();
    }

    inline double f(const vect3& x) const
    {
        //Analytical value of the inner integral in operator D. See DeMunck article for further details.
        //  for non-optimized version of operator D
        //  returns the value of the inner integral of operator D on a triangle used for a P1 function
        vect3 Y1=v1-x;
        vect3 Y2=v2-x;
        vect3 Y3=v3-x;
        double y1=Y1.norme();
        double y2=Y2.norme();
        double y3=Y3.norme();
        double d=Y1*(Y2^Y3);

        double derr=1e-10;
        if(fabs(d)<derr) return 0.0;

        double omega=2*atan2(d,(y1*y2*y3+y1*(Y2*Y3)+y2*(Y3*Y1)+y3*(Y1*Y2)));

        vect3 Z1=Y2^Y3;
        vect3 Z2=Y3^Y1;
        vect3 Z3=Y1^Y2;
        vect3 D1=Y2-Y1;
        vect3 D2=Y3-Y2;
        vect3 D3=Y1-Y3;
        double d1=D1.norme();
        double d2=D2.norme();
        double d3=D3.norme();
        double g1=-1/d1*log((y1*d1+Y1*D1)/(y2*d1+Y2*D1));
        double g2=-1/d2*log((y2*d2+Y2*D2)/(y3*d2+Y3*D2));
        double g3=-1/d3*log((y3*d3+Y3*D3)/(y1*d3+Y1*D3));
        vect3 N=Z1+Z2+Z3;
        double A=N.norme2();
        vect3 S=D1*g1+D2*g2+D3*g3;
        double omega_i[3];
        omega_i[0]=1/A*(Z1*N*omega+d*(D2*S));
        omega_i[1]=1/A*(Z2*N*omega+d*(D3*S));
        omega_i[2]=1/A*(Z3*N*omega+d*(D1*S));
        double result=omega_i[i-1];

        return result;
    }
};


class analytiqueD3 : public fContainer<vect3> {
private:
    vect3 v1,v2,v3;
    double aire;
public:
    analytiqueD3(){}
    ~analytiqueD3(){}
    inline void init( const mesh& m1, const int trg)
    {
        v1=m1.vctr(m1.trngl(trg).id1());
        v2=m1.vctr(m1.trngl(trg).id2());
        v3=m1.vctr(m1.trngl(trg).id3());
        aire=m1.trngl(trg).getAire();
    }

    inline vect3 f(const vect3& x) const
    {
        //Analytical value of the inner integral in operator D. See DeMunck article for further details.
        //  for non-optimized version of operator D
        //  returns in a vector, the inner integrals of operator D on a triangle viewed as a part of the 3
        //  P1 functions it has a part in.
        vect3 Y1=v1-x;
        vect3 Y2=v2-x;
        vect3 Y3=v3-x;
        double y1=Y1.norme();
        double y2=Y2.norme();
        double y3=Y3.norme();
        double d=Y1*(Y2^Y3);

        double derr=1e-10;
        if(fabs(d)<derr) return 0.0;

        double omega=2*atan2(d,(y1*y2*y3+y1*(Y2*Y3)+y2*(Y3*Y1)+y3*(Y1*Y2)));

        vect3 Z1=Y2^Y3;
        vect3 Z2=Y3^Y1;
        vect3 Z3=Y1^Y2;
        vect3 D1=Y2-Y1;
        vect3 D2=Y3-Y2;
        vect3 D3=Y1-Y3;
        double d1=D1.norme();
        double d2=D2.norme();
        double d3=D3.norme();
        double g1=-1/d1*log((y1*d1+Y1*D1)/(y2*d1+Y2*D1));
        double g2=-1/d2*log((y2*d2+Y2*D2)/(y3*d2+Y3*D2));
        double g3=-1/d3*log((y3*d3+Y3*D3)/(y1*d3+Y1*D3));
        vect3 N=Z1+Z2+Z3;
        double A=N.norme2();
        vect3 S=D1*g1+D2*g2+D3*g3;
        vect3 omega_i;
        omega_i._x()=1/A*(Z1*N*omega+d*(D2*S));
        omega_i._y()=1/A*(Z2*N*omega+d*(D3*S));
        omega_i._z()=1/A*(Z3*N*omega+d*(D1*S));

        return omega_i;
    }
};

class analytiqueDipPot: public fContainer<double> {
private:
    vect3 q,r0;
public:
    analytiqueDipPot(){}
    ~analytiqueDipPot(){}

    inline void init( const vect3 &_q, const vect3 &_r0)
    {
        q=_q;
        r0=_r0;
    }

    inline double f(const vect3& x) const
    {
        vect3 r0r=x-r0;
        double r0rn3=r0r.norme(); r0rn3=r0rn3*r0rn3*r0rn3;
        return (q*r0r)/r0rn3;
    }
};

class analytiqueDipPotDer : public fContainer<vect3> {
private:
    vect3 q,r0;
    vect3 H0,H1,H2;
    vect3 H0p0DivNorm2,H1p1DivNorm2,H2p2DivNorm2,n;
public:
    analytiqueDipPotDer(){}
    ~analytiqueDipPotDer(){}
    inline void init( const mesh& m, const int nT, const vect3 &_q, const vect3 _r0)
    {
        q=_q;
        r0=_r0;

        triangle& T=(triangle&)m.trngl(nT);

        vect3 p0,p1,p2,p1p0,p2p1,p0p2,p1p0n,p2p1n,p0p2n,p1H0,p2H1,p0H2;
        p0=m.vctr(T.som1());
        p1=m.vctr(T.som2());
        p2=m.vctr(T.som3());


        p1p0=p0-p1; p2p1=p1-p2; p0p2=p2-p0;
        p1p0n=p1p0; p1p0n.normalize(); p2p1n=p2p1; p2p1n.normalize(); p0p2n=p0p2; p0p2n.normalize();

        p1H0=(p1p0*p2p1n)*p2p1n; H0=p1H0+p1; H0p0DivNorm2=p0-H0; H0p0DivNorm2=H0p0DivNorm2/H0p0DivNorm2.norme2();
        p2H1=(p2p1*p0p2n)*p0p2n; H1=p2H1+p2; H1p1DivNorm2=p1-H1; H1p1DivNorm2=H1p1DivNorm2/H1p1DivNorm2.norme2();
        p0H2=(p0p2*p1p0n)*p1p0n; H2=p0H2+p0; H2p2DivNorm2=p2-H2; H2p2DivNorm2=H2p2DivNorm2/H2p2DivNorm2.norme2();

        n=-p1p0^p0p2;
        n.normalize();

    }

    inline vect3 f(const vect3& x) const
    {
        vect3 rr0=r0-x;
        vect3 P1part(H0p0DivNorm2*(x-H0),H1p1DivNorm2*(x-H1),H2p2DivNorm2*(x-H2));
        double inv_r0rn3=(rr0).norme(); inv_r0rn3=inv_r0rn3*inv_r0rn3*inv_r0rn3; inv_r0rn3=1.0/inv_r0rn3;
        rr0.normalize();
        double EMpart=q*n*inv_r0rn3-3*inv_r0rn3*(q*rr0)*(n*rr0);
        return EMpart*P1part;
    }
};

#endif
