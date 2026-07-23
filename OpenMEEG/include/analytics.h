// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <cmath>
#include <triangle.h>
#include <dipole.h>

namespace OpenMEEG {

    // Analytical computation of the integral over a triangle of 1/|r-x| dr.
    // The formula are derived from Ferguson, Zhang and Stroink, 1994.

    class OPENMEEG_EXPORT OperatorS {
    public:

        OperatorS(const Triangle& T): triangle(T) {
            for (unsigned i=0; i<3; ++i) {
                const Vect3  vref = triangle.vertex(i)-triangle.vertex(next[i]);
                const double len  = vref.norm();
                const Vect3  unit = vref/len;
                lengths[i] = len;
                normals[i] = crossprod(triangle.normal(),unit);
            }
        }

        // Analytical value of the internal integral of S operator at point x.

        double operator()(const Vect3& x) const {
            const Vect3  rays[3]  = { triangle.vertex(0)-x, triangle.vertex(1)-x, triangle.vertex(2)-x };
            const double dists[3] = { rays[0].norm(),       rays[1].norm(),       rays[2].norm()       };
            const double R[3]     = { dists[0]+dists[1],    dists[1]+dists[2],    dists[2]+dists[0]    };

            // Formula derived from the article Ferguson, Zhang, Stroink IEEE, Trans. on Biomedical Engineering (41) 5, 1994.
            // The result in the article is (n,p0-x,p1-x)/|p0-p1| log(arg) with:
            // arg = (|p1-x|*|p1-p0|+dotprod(p1-x,p1-p0))/(|p0-x|*|p1-p0|+dotprod(p0-x,p1-p0));
            // It is not too difficult (but a bit tedious) to prove that:
            // arg = (|p0-x|+|p1-x|+|p1-p0|)/(|p0-x|+|p1-x|-|p1-p0|);
            // This last form of arg is both less expensive and makes it clear that arg>=1,
            // The lengths |p0-x|,|p1-x| and |p1-p0| are respectively dists[0], dists[1] and lengths[0].
            // |p0-x|+|p1-x| is R[0]. Introducing u = |p1-p0|/(|p0-x|+|p1-x|) leads to:
            // arg = (1+u)/(1-u) and log(arg) = 2 atanh(u)
            // This last formula is both more stable and simpler to compute.

            double sum = 0.0;
            for (unsigned i=0; i<3; ++i)
                sum += dotprod(rays[i],normals[i])*atanh(lengths[i]/R[i]);

            return 2*sum-dotprod(rays[0],triangle.normal())*triangle.solid_angle(x);
        }

    private:

        static constexpr unsigned next[3] = { 1, 2, 0 };

        const Triangle& triangle;

        // Storage order in lengths and normals is p1p0, p2p1 and p0p2 (p[next[i]]-p[i], for i=0,1,2, where p are the vertices of the triangle).

        double lengths[3]; ///< Lengths of the edges of the triangle. 
        Vect3  normals[3]; ///< These are the normals in the triangle plane to the edges of the triangle.
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
