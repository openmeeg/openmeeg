/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

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

#pragma once

#include <isnormal.H>
#include <triangle.h>
#include <dipole.h>

namespace OpenMEEG {

    inline double integral_simplified_green(const Vect3& p0x, const double norm2p0x,
                                            const Vect3& p1x, const double norm2p1x,
                                            const Vect3& p1p0,const double norm2p1p0)
    {
        //  The quantity arg is normally >= 1, verifying this relates to a triangular inequality
        //  between p0, p1 and x.
        //  Consequently, there is no need of an absolute value in the first case.

        const double arg = (norm2p0x*norm2p1p0-dotprod(p0x,p1p0))/(norm2p1x*norm2p1p0-dotprod(p1x,p1p0));
        return (std::isnormal(arg) && arg>0.0) ? log(arg) : fabs(log(norm2p1x/norm2p0x));
    }

    class OPENMEEG_EXPORT analyticS {

        void initialize(const Vect3& v0,const Vect3& v1,const Vect3& v2) {
            // All computations needed when the first triangle of integration is changed

            p0 = v0;
            p1 = v1;
            p2 = v2;

            p1p0 = p1-p0;
            p2p1 = p2-p1;
            p0p2 = p0-p2;

            norm2p1p0 = p1p0.norm();
            norm2p2p1 = p2p1.norm();
            norm2p0p2 = p0p2.norm();
        }

        void finish_intialization() {
            nu0 = (p1p0^n);
            nu1 = (p2p1^n);
            nu2 = (p0p2^n);
            nu0.normalize();
            nu1.normalize();
            nu2.normalize();
        }

    public:

        analyticS(const Triangle& T) {
            initialize(T.vertex(0),T.vertex(1),T.vertex(2));
            n = T.normal();
            finish_intialization();
        }

        analyticS(const Vect3& v0,const Vect3& v1,const Vect3& v2) {
            initialize(v0,v1,v2);
            n = p1p0^p0p2;
            n /= n.norm();
            finish_intialization();
        }

        double f(const Vect3& x) const {
            // analytical value of the internal integral of S operator at point X
            const Vect3& p0x = p0-x;
            const Vect3& p1x = p1-x;
            const Vect3& p2x = p2-x;
            const double norm2p0x = p0x.norm();
            const double norm2p1x = p1x.norm();
            const double norm2p2x = p2x.norm();

            const double g0 = integral_simplified_green(p0x,norm2p0x,p1x,norm2p1x,p1p0,norm2p1p0);
            const double g1 = integral_simplified_green(p1x,norm2p1x,p2x,norm2p2x,p2p1,norm2p2p1);
            const double g2 = integral_simplified_green(p2x,norm2p2x,p0x,norm2p0x,p0p2,norm2p0p2);

            const double alpha = dotprod(p0x,n);

            return ((dotprod(p0x,nu0)*g0+dotprod(p1x,nu1)*g1+dotprod(p2x,nu2)*g2)-alpha*x.solid_angle(p0,p1,p2));
        }

    private:

        Vect3 p0, p1, p2; //!< vertices of the triangle
        Vect3 p2p1, p1p0, p0p2;
        Vect3 nu0, nu1, nu2;
        Vect3 n;
        double norm2p2p1, norm2p1p0, norm2p0p2;
        double tanTHETA0m, tanTHETA0p, tanTHETA1m, tanTHETA1p, tanTHETA2m, tanTHETA2p;
    };

    class OPENMEEG_EXPORT analyticD3 {

        static Vect3 unit_vector(const Vect3& V) { return V/V.norm(); }

        Vect3 diff(const unsigned i,const unsigned j) const { return triangle.vertex(i)-triangle.vertex(j); }

        // TODO: Introduce a D matrix, and a dipole....

        #if 0
        Matrix initD() const {
            Matrix res(3,3);
            res.setlin(1,D2);
            res.setlin(2,D3);
            res.setlin(3,D1);
            return res;
        }
        #endif

    public:

        analyticD3(const Triangle& T):
            triangle(T),D1(diff(1,0)),D2(diff(2,1)),D3(diff(0,2)),U1(unit_vector(D1)),U2(unit_vector(D2)),U3(unit_vector(D3))
        { }

        inline Vect3 f(const Vect3& x) const {
            // Analytical value of the inner integral in operator D. See DeMunck article for further details.
            // Used in non-optimized version of operator D.
            // Returns a vector of the inner integrals of operator D on a triangle wrt its three P1 functions.

            //  First part omega is just x.solid_angle(v1,v2,v3)

            const Vect3& Y1 = triangle.vertex(0)-x;
            const Vect3& Y2 = triangle.vertex(1)-x;
            const Vect3& Y3 = triangle.vertex(2)-x;
            const double y1 = Y1.norm();
            const double y2 = Y2.norm();
            const double y3 = Y3.norm();
            const double d  = det(Y1,Y2,Y3);

            if (fabs(d)<1e-10)
                return 0.0;

            const double omega = 2*atan2(d,(y1*y2*y3+y1*dotprod(Y2,Y3)+y2*dotprod(Y3,Y1)+y3*dotprod(Y1,Y2)));

            const Vect3& Z1 = crossprod(Y2,Y3);
            const Vect3& Z2 = crossprod(Y3,Y1);
            const Vect3& Z3 = crossprod(Y1,Y2);
            const double g1 = log((y2+dotprod(Y2,U1))/(y1+dotprod(Y1,U1)));
            const double g2 = log((y3+dotprod(Y3,U2))/(y2+dotprod(Y2,U2)));
            const double g3 = log((y1+dotprod(Y1,U3))/(y3+dotprod(Y3,U3)));
            const Vect3& N = Z1+Z2+Z3;
            const double invA = 1.0/N.norm2();
            const Vect3& S = U1*g1+U2*g2+U3*g3;

            return (omega*Vect3(dotprod(Z1,N),dotprod(Z2,N),dotprod(Z3,N))+d*Vect3(dotprod(D2,S),dotprod(D3,S),dotprod(D1,S)))/N.norm2();
        }

    private:

        const Triangle& triangle;
        //const Matrix    D;
        const Vect3     D1;
        const Vect3     D2;
        const Vect3     D3;
        const Vect3     U1;
        const Vect3     U2;
        const Vect3     U3;
    };

    class OPENMEEG_EXPORT analyticDipPotDer {
    public:

        analyticDipPotDer(const Dipole& dip,const Triangle& T): dipole(dip) {

            const Vect3& p0 = T.vertex(0);
            const Vect3& p1 = T.vertex(1);
            const Vect3& p2 = T.vertex(2);

            const Vect3& p1p0 = p0-p1;
            const Vect3& p2p1 = p1-p2;
            const Vect3& p0p2 = p2-p0;
            const Vect3& p1p0n = p1p0/p1p0.norm();
            const Vect3& p2p1n = p2p1/p2p1.norm();
            const Vect3& p0p2n = p0p2/p0p2.norm();

            const Vect3& p1H0 = dotprod(p1p0,p2p1n)*p2p1n;
            H0 = p1H0+p1;
            H0p0DivNorm2 = p0-H0;
            H0p0DivNorm2 = H0p0DivNorm2/H0p0DivNorm2.norm2();
            const Vect3& p2H1 = dotprod(p2p1,p0p2n)*p0p2n;
            H1 = p2H1+p2;
            H1p1DivNorm2 = p1-H1;
            H1p1DivNorm2 = H1p1DivNorm2/H1p1DivNorm2.norm2();
            const Vect3& p0H2 = dotprod(p0p2,p1p0n)*p1p0n;
            H2 = p0H2+p0;
            H2p2DivNorm2 = p2-H2;
            H2p2DivNorm2 = H2p2DivNorm2/H2p2DivNorm2.norm2();

            n = -crossprod(p1p0,p0p2);
            n.normalize();
        }

        Vect3 f(const Vect3& r) const {
            Vect3 P1part(dotprod(H0p0DivNorm2,r-H0),dotprod(H1p1DivNorm2,r-H1),dotprod(H2p2DivNorm2,r-H2));

            // B = n.grad_x(A) with grad_x(A)= q/||^3 - 3r(q.r)/||^5

            const Vect3& x         = r-dipole.position();
            const double inv_xnrm2 = 1.0/x.norm2();
            const double EMpart = dotprod(n,dipole.moment()-3*dotprod(dipole.moment(),x)*x*inv_xnrm2)*(inv_xnrm2*sqrt(inv_xnrm2));

            return -EMpart*P1part; // RK: why - sign ?
        }

    private:

        const Dipole& dipole;

        Vect3 H0, H1, H2;
        Vect3 H0p0DivNorm2, H1p1DivNorm2, H2p2DivNorm2, n;
    };
}
