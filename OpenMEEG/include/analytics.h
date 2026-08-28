// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <cmath>
#include <triangle.h>
#include <BarycentricCoordinates.h>

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

        double operator()(const Point3D& x) const {
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

    // Analytical computation of the integral over a triangle T of the P1 single-layer potential (integral of lambda_i(y)/|y-x| over y).
    // The implementation follows a formula obtained in Graglia 1993.
    // Used in non-optimized version of operator D.
    // Returns a vector of the inner integrals of operator D on a triangle w.r.t. its three P1 functions.
    // See DeMunck article for further details.

    class OPENMEEG_EXPORT P1SingleLayerPotential {

        Vect3 diff(const unsigned i,const unsigned j) const { return triangle.vertex(i)-triangle.vertex(j); }

    public:

        P1SingleLayerPotential(const Triangle& T):
            triangle(T),
            D{ diff(1,0), diff(2,1), diff(0,2)},
            U{ D[0].unit_vector(), D[1].unit_vector(), D[2].unit_vector() }
        { }

        Vect3 operator()(const Point3D& x) const {

            const Vect3 rays[3]  = { triangle.vertex(0)-x,       triangle.vertex(1)-x,       triangle.vertex(2)-x       };
            const Vect3 dists    = { rays[0].norm(),             rays[1].norm(),             rays[2].norm()             };
            const Vect3 Z[3]     = { crossprod(rays[1],rays[2]), crossprod(rays[2],rays[0]), crossprod(rays[0],rays[1]) };
            const Vect3 prods    = { dotprod(rays[1],rays[2]),   dotprod(rays[2],rays[0]),   dotprod(rays[0],rays[1])   };

            const Vect3& n = triangle.normal();
            const double h = dotprod(rays[0],triangle.normal()); // Half the distance of x to the plane of the triangle.
            const double d = 2*h*triangle.area();

            // If the volume of the tetrahedron (x,p0,p1,p2), i.e. h is too small just return 0.

            if (fabs(h)<1e-10)
                return 0.0;

            Vect3 S;
            for (unsigned i=0; i<3; ++i) {
                const unsigned i1 = next[i];
                S  += U[i]*log((dists[i1]+dotprod(rays[i1],U[i]))/(dists[i]+dotprod(rays[i],U[i])));
            }

            // omega is just half of triangle.solid_angle(x).
            // We duplicated the code here because many quantities are needed for the rest of the computation.

            const double d1    = dists[0]*dists[1]*dists[2]+dotprod(dists,prods);
            const double omega = atan2(d,d1);

            const Vect3 E1(dotprod(Z[0],n),dotprod(Z[1],n),dotprod(Z[2],n));
            const Vect3 E2(dotprod(D[1],S),dotprod(D[2],S),dotprod(D[0],S));

            return (omega*E1+0.5*h*E2)/triangle.area();
        }

    private:

        static constexpr unsigned next[3] = { 1, 2, 0 };

        const Triangle& triangle;
        const Vect3     D[3];
        const Vect3     U[3];
    };

    // P1 integral of the normal component of the gradient of the potential.

    template <typename Source>
    class OPENMEEG_EXPORT PotentialDerivative {
    public:

        PotentialDerivative(const Source& src,const Triangle& T):
            barycentric_coords(T),triangle(T),source(src)
        { }

        Vect3 operator()(const Point3D& r) const {

            // In this case, r is a point in triangle.
            // We compute the barycentric coordinates of r in the triangle.
            // These are exactly the P1 coordinates of the point in the triangle.

            const double pot_grad_n = source.gradient(r,triangle.normal());

            return pot_grad_n*barycentric_coords(r);
        }

    private:

        const BarycentricCoordinates barycentric_coords;

        const Triangle& triangle;
        const Source&   source;
    };
}
