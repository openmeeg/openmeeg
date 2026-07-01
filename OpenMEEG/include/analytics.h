// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <cmath>
#include <triangle.h>
#include <monopole.h>
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

    // TODO: There is some common code in analyticS and analyticD3 (compute distances of the point to the triangle vertices),
    // even though it takes a somewhat different form.
    // Factorize this ?

    class OPENMEEG_EXPORT analyticS {

        // TODO: Instead of storing the individual points and the normal, store the triangle instead.

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

            //  First part omega is just x.solid_angle(triangle.vertex(0),triangle.vertex(1),triangle.vertex(2))

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

    // This class basically recomputes the barycentric coordinates of a point r in the triangle.
    // This just solves with a least squares approach the linear system
    //                 r = l*r0+m*r1+(1-l-m)*r2
    // Least squares is needed because if we select a subset of equations, we may run in a
    // singular case. This leads to the 2x2 system:
    //
    //              |r0-r2]^2       * l + (r0-r2).(r1-r1) * m = (r0-r2).(r-r2)
    //              (r0-r2).(r1-r1) * l + |r1-r2]^2       * m = (r1-r2).(r-r2)
    //
    //  which is solved using the Kramer formualas.             
    //
    // TODO: It is a little bit silly to do this. The points that are passed by the integrator
    // are computed from their barycentric coordinates and we now do the opposite computation.
    // A more clever way would be to pass the barycentric coordinates directly and to compute
    // the point in the analytic functions. This will totally suppress the need for this class.

    class OPENMEEG_EXPORT BarycentricCoordinates {
    public:

        BarycentricCoordinates(const Triangle& T) {

            const Vect3& p0 = T.vertex(0);
            const Vect3& p1 = T.vertex(1);
            const Vect3& p2 = T.vertex(2);

            r20 = p0-p2;
            r21 = p1-p2;
            r2  = p2;

            const double norm_r20_squared = r20.norm2();
            const double norm_r21_squared = r21.norm2();
            const double prod_r20_r21     = dotprod(r20,r21);

            const double det = norm_r20_squared*norm_r21_squared-sqr(prod_r20_r21);
            const double invdet = 1.0/det;

            c11 = norm_r20_squared*invdet;
            c22 = norm_r21_squared*invdet;
            c12 = prod_r20_r21*invdet;
        }

        Vect3 operator()(const Vect3& r) const {

            const Vect3& u = r-r2;
            const double b1 = dotprod(r20,u);
            const double b2 = dotprod(r21,u);
            const double l = b1*c22-b2*c12;
            const double m = b2*c11-b1*c12;
            return Vect3(l,m,1-l-m);
        }

    private:

        double c11,c12,c22;
        Vect3  r20,r21,r2;
    };

    class OPENMEEG_EXPORT analyticMonopolePotDer {
    public:

        analyticMonopolePotDer(const Monopole& monop,const Triangle& T):
            barycentric_coords(T),triangle(T),monopole(monop)
        { }

        Vect3 f(const Vect3& r) const {

            // B = n.grad_x(A) with grad_x(A)= -q x/|x|^3, with x = r-r0

            const Vect3& x      = r-monopole.position();
            const double xnrm2  = x.norm2();
            const Vect3& n      = triangle.normal();
			const double EMpart = monopole.charge()*dotprod(n,x)/(xnrm2*sqrt(xnrm2));
            const Vect3& P1term = barycentric_coords(r); // P1 function values are just the barycentric coordinates of r.

            return EMpart*P1term; // RK: why not - sign ?
        }

    private:

        const BarycentricCoordinates barycentric_coords;

        const Triangle& triangle;
        const Monopole& monopole;
    };

    class OPENMEEG_EXPORT analyticDipPotDer {
    public:

        analyticDipPotDer(const Dipole& dip,const Triangle& T):
            barycentric_coords(T),triangle(T),dipole(dip)
        { }

        Vect3 f(const Vect3& r) const {

            // B = n.grad_x(A) with grad_x(A)= q/|x|^3 - 3r(q.x)/|x|^5 with x = r-r0

            const Vect3& x         = r-dipole.position();
            const double inv_xnrm2 = 1.0/x.norm2();
            const Vect3& n         = triangle.normal();
            const double EMpart    = dotprod(n,dipole.moment()-(3*dotprod(dipole.moment(),x)*inv_xnrm2)*x)*(inv_xnrm2*sqrt(inv_xnrm2));
            const Vect3& P1term    = barycentric_coords(r); // P1 function values are just the barycentric coordinates of r.

            return -EMpart*P1term; // RK: why - sign ?
        }

    private:

        const BarycentricCoordinates barycentric_coords;

        const Triangle& triangle;
        const Dipole&   dipole;
    };
}
