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

#ifndef OPENMEEG_ANALYTICS_H
#define OPENMEEG_ANALYTICS_H

#include <stdio.h>

#include "mesh3.h"

namespace OpenMEEG {

    class OPENMEEG_EXPORT analyticS {
    private:
        Vect3 p0,p1,p2; //!< vertices of the triangle
        Vect3 p2p1,p1p0,p0p2;
        Vect3 nu0,nu1,nu2;
        Vect3 n;
        double norm2p2p1,norm2p1p0,norm2p0p2;
        double tanTHETA0m,tanTHETA0p,tanTHETA1m,tanTHETA1p,tanTHETA2m,tanTHETA2p;
        
        void init_aux() {
            n *= -(1/n.norm()) ;
            nu0 = (n^p1p0);
            nu1 = (n^p2p1);
            nu2 = (n^p0p2);
            nu0.normalize();
            nu1.normalize();
            nu2.normalize();
        }

    public:

        analyticS(){}
        ~analyticS(){}
        void init( int nT, const Mesh &m )
        {
            // all computations needed when the first triangle of integration is changed
            Triangle& T = (Triangle&) m.getTrg(nT);

            p0 = m.getPt(T.s1());
            p1 = m.getPt(T.s2());
            p2 = m.getPt(T.s3());

            p1p0 = p1-p0; p2p1 = p2-p1; p0p2 = p0-p2;
            norm2p1p0 = p1p0.norm(); norm2p2p1 = p2p1.norm(); norm2p0p2 = p0p2.norm();

            n = T.normal();
            init_aux();
        }

        void init( const Vect3& v0, const Vect3& v1, const Vect3& v2 )
        {
            // all computations needed when the first triangle of integration is changed
            p0 = v0;
            p1 = v1;
            p2 = v2;
            p1p0 = p1-p0; p2p1 = p2-p1; p0p2 = p0-p2;
            norm2p1p0 = p1p0.norm(); norm2p2p1 = p2p1.norm(); norm2p0p2 = p0p2.norm();

            n = p1p0^p0p2;
            init_aux();
        }

        inline double f(const Vect3& x) const
        {
            // analytical value of the internal integral of S operator at point X
            Vect3 p1x = p1-x, p2x = p2-x, p0x = p0-x ;
            double norm2p0x = p0x.norm();
            double norm2p1x = p1x.norm();
            double norm2p2x = p2x.norm();
            double alpha = (x-p0)*n ;
            double g0,g1,g2;

            if ((p0x^p1p0).norm() > .00000001)
                g0 = -log(norm2p1x-p1x*p1p0*(1.0/norm2p1p0) )+log(norm2p0x-p0x*p1p0*(1.0/norm2p1p0) );
            else
                g0= fabs(log(norm2p1x)-log(norm2p0x));

            if ((p1x^p2p1).norm() > .00000001)
                g1 = -log(norm2p2x-p2x*p2p1*(1.0/norm2p2p1) )+log(norm2p1x-p1x*p2p1*(1.0/norm2p2p1) );
            else
                g1= fabs(log(norm2p2x)-log(norm2p1x));

            if ((p2x^p0p2).norm() > .00000001)
                g2 = -log(norm2p0x-p0x*p0p2*(1.0/norm2p0p2) )+log(norm2p2x-p2x*p0p2*(1.0/norm2p0p2) );
            else
                g2 = fabs(log(norm2p0x)-log(norm2p2x));

            return ((p0x*nu0)*g0+(p1x*nu1)*g1+(p2x*nu2)*g2)-alpha*x.solangl(p0,p1,p2);
        }
    };

    class OPENMEEG_EXPORT analyticD {
    private:
        Vect3 v1,v2,v3;
        int i;
        double aire;
    public:
        analyticD(){}
        ~analyticD(){}
        inline void init( const Mesh& m1, const int trg, const int noeud)
        {
            v1 = m1.getPt(m1.getTrg(trg).s1());
            v2 = m1.getPt(m1.getTrg(trg).s2());
            v3 = m1.getPt(m1.getTrg(trg).s3());
            i = m1.getTrg(trg).contains(noeud);
            aire = m1.getTrg(trg).getArea();
        }

        inline double f(const Vect3& x) const
        {
            //Analytical value of the inner integral in operator D. See DeMunck article for further details.
            //  for non-optimized version of operator D
            //  returns the value of the inner integral of operator D on a triangle used for a P1 function
            Vect3 Y1 = v1-x;
            Vect3 Y2 = v2-x;
            Vect3 Y3 = v3-x;
            double y1 = Y1.norm();
            double y2 = Y2.norm();
            double y3 = Y3.norm();
            double d = Y1*(Y2^Y3);

            double derr = 1e-10;
            if(fabs(d)<derr) return 0.0;

            double omega = 2*atan2(d,(y1*y2*y3+y1*(Y2*Y3)+y2*(Y3*Y1)+y3*(Y1*Y2)));

            Vect3 Z1 = Y2^Y3;
            Vect3 Z2 = Y3^Y1;
            Vect3 Z3 = Y1^Y2;
            Vect3 D1 = Y2-Y1;
            Vect3 D2 = Y3-Y2;
            Vect3 D3 = Y1-Y3;
            double d1 = D1.norm();
            double d2 = D2.norm();
            double d3 = D3.norm();
            double g1 = -1/d1*log((y1*d1+Y1*D1)/(y2*d1+Y2*D1));
            double g2 = -1/d2*log((y2*d2+Y2*D2)/(y3*d2+Y3*D2));
            double g3 = -1/d3*log((y3*d3+Y3*D3)/(y1*d3+Y1*D3));
            Vect3 N = Z1+Z2+Z3;
            double A = N.norm2();
            Vect3 S = D1*g1+D2*g2+D3*g3;
            double omega_i[3];
            omega_i[0] = 1/A*(Z1*N*omega+d*(D2*S));
            omega_i[1] = 1/A*(Z2*N*omega+d*(D3*S));
            omega_i[2] = 1/A*(Z3*N*omega+d*(D1*S));
            double result = omega_i[i-1];

            return result;
        }
    };


    class OPENMEEG_EXPORT analyticD3 {
    private:
        Vect3 v1,v2,v3;
        double aire;
    public:
        analyticD3(){}
        ~analyticD3(){}
        inline void init( const Mesh& m1, const int trg)
        {
            v1 = m1.getPt(m1.getTrg(trg).s1());
            v2 = m1.getPt(m1.getTrg(trg).s2());
            v3 = m1.getPt(m1.getTrg(trg).s3());
            aire = m1.getTrg(trg).getArea();
        }

        inline Vect3 f(const Vect3& x) const
        {
            //Analytical value of the inner integral in operator D. See DeMunck article for further details.
            //  for non-optimized version of operator D
            //  returns in a vector, the inner integrals of operator D on a triangle viewed as a part of the 3
            //  P1 functions it has a part in.
            Vect3 Y1 = v1-x;
            Vect3 Y2 = v2-x;
            Vect3 Y3 = v3-x;
            double y1 = Y1.norm();
            double y2 = Y2.norm();
            double y3 = Y3.norm();
            double d = Y1*(Y2^Y3);

            double derr = 1e-10;
            if(fabs(d)<derr) return 0.0;

            double omega = 2*atan2(d,(y1*y2*y3+y1*(Y2*Y3)+y2*(Y3*Y1)+y3*(Y1*Y2)));

            Vect3 Z1 = Y2^Y3;
            Vect3 Z2 = Y3^Y1;
            Vect3 Z3 = Y1^Y2;
            Vect3 D1 = Y2-Y1;
            Vect3 D2 = Y3-Y2;
            Vect3 D3 = Y1-Y3;
            double d1 = D1.norm();
            double d2 = D2.norm();
            double d3 = D3.norm();
            double g1 = -1/d1*log((y1*d1+Y1*D1)/(y2*d1+Y2*D1));
            double g2 = -1/d2*log((y2*d2+Y2*D2)/(y3*d2+Y3*D2));
            double g3 = -1/d3*log((y3*d3+Y3*D3)/(y1*d3+Y1*D3));
            Vect3 N = Z1+Z2+Z3;
            double A = N.norm2();
            Vect3 S = D1*g1+D2*g2+D3*g3;
            Vect3 omega_i;
            omega_i.x() = 1/A*(Z1*N*omega+d*(D2*S));
            omega_i.y() = 1/A*(Z2*N*omega+d*(D3*S));
            omega_i.z() = 1/A*(Z3*N*omega+d*(D1*S));

            return omega_i;
        }
    };

    class OPENMEEG_EXPORT analyticDipPot {
    private:
        Vect3 q,r0;
    public:
        analyticDipPot(){}
        ~analyticDipPot(){}

        inline void init( const Vect3& _q, const Vect3& _r0)
        {
            q = _q;
            r0 = _r0;
        }

        inline double f(const Vect3& x) const
        {
            // RK: A = q.(x-r0)/||^3
            Vect3 r = x-r0;
            double rn = r.norm();
            return (q*r)/(rn*rn*rn);
        }
    };

    class OPENMEEG_EXPORT analyticDipPotDer {
    private:
        Vect3 q,r0;
        Vect3 H0,H1,H2;
        Vect3 H0p0DivNorm2,H1p1DivNorm2,H2p2DivNorm2,n;
    public:
        analyticDipPotDer(){}
        ~analyticDipPotDer(){}
        inline void init( const Mesh& m, const int nT, const Vect3 &_q, const Vect3 _r0)
        {
            q = _q;
            r0 = _r0;

            Triangle& T = (Triangle&)m.getTrg(nT);

            Vect3 p0,p1,p2,p1p0,p2p1,p0p2,p1p0n,p2p1n,p0p2n,p1H0,p2H1,p0H2;
            p0 = m.getPt(T.s1());
            p1 = m.getPt(T.s2());
            p2 = m.getPt(T.s3());

            p1p0 = p0-p1; p2p1 = p1-p2; p0p2 = p2-p0;
            p1p0n = p1p0; p1p0n.normalize(); p2p1n = p2p1; p2p1n.normalize(); p0p2n = p0p2; p0p2n.normalize();

            p1H0 = (p1p0*p2p1n)*p2p1n; H0 = p1H0+p1; H0p0DivNorm2 = p0-H0; H0p0DivNorm2 = H0p0DivNorm2/H0p0DivNorm2.norm2();
            p2H1 = (p2p1*p0p2n)*p0p2n; H1 = p2H1+p2; H1p1DivNorm2 = p1-H1; H1p1DivNorm2 = H1p1DivNorm2/H1p1DivNorm2.norm2();
            p0H2 = (p0p2*p1p0n)*p1p0n; H2 = p0H2+p0; H2p2DivNorm2 = p2-H2; H2p2DivNorm2 = H2p2DivNorm2/H2p2DivNorm2.norm2();

            n = -p1p0^p0p2;
            n.normalize();
        }

        inline Vect3 f(const Vect3& x) const
        {
            Vect3 P1part(H0p0DivNorm2*(x-H0),H1p1DivNorm2*(x-H1),H2p2DivNorm2*(x-H2));

            // RK: B = n.grad_x(A) with grad_x(A)= q/||^3 - 3r(q.r)/||^5
            Vect3 r = x-r0;
            double rn = r.norm();
            double EMpart = n*(q/pow(rn,3.)-3*(q*r)*r/pow(rn,5.));

            return -EMpart*P1part; // RK: why - sign ?
        }
    };
}

#endif  //! OPENMEEG_ANALYTICS_H
