// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <triangle.h>

namespace OpenMEEG {

    // This class basically recomputes the barycentric coordinates of a point r in the triangle.
    // This just solves with a least squares approach the linear system
    //                 r = l*r0+m*r1+(1-l-m)*r2
    // Least squares is needed because if we select a subset of equations, we may run in a
    // singular case. This leads to the 2x2 system:
    //
    //              |r0-r2]^2       * l + (r0-r2).(r1-r1) * m = (r0-r2).(r-r2)
    //              (r0-r2).(r1-r1) * l + |r1-r2]^2       * m = (r1-r2).(r-r2)
    //
    //  which is solved using the Kramer formulas.             
    //
    // TODO: It is a little bit silly to do this. The points that are passed by the integrator
    // are computed from their barycentric coordinates and we now do the opposite computation.
    // A more clever way would be to pass the barycentric coordinates directly and to compute
    // the point in the analytic functions. This will totally suppress the need for this class.

    class OPENMEEG_EXPORT BarycentricCoordinates {
    public:

        BarycentricCoordinates(const Triangle& T): triangle(T) {

            edge10 = triangle.vertex(1)-triangle.vertex(0);
            edge20 = triangle.vertex(2)-triangle.vertex(0);

            // Coefficients of the 2x2 symmetric matrix, for finding the barycentric coordinates.

            const double c11 = dotprod(edge10,edge10);
            const double c12 = dotprod(edge10,edge20);
            const double c22 = dotprod(edge20,edge20);

            const double inv_det = 1.0/(c11*c22-c12*c12);

            // Coefficients of the inverse symmetric matrix.

            d11 =  c22*inv_det;
            d12 = -c12*inv_det;
            d22 =  c11*inv_det;
        }

        Vect3 operator()(const Vect3& r) const {

            // Barycentric coordinates of r in the triangle.

            const Vect3& v   = r-triangle.vertex(0);
            const double d1 = dotprod(v,edge10);
            const double d2 = dotprod(v,edge20);

            // Apply the matrix [[d11,d12],[d12,d22]] to vector [d1,d2]^T to obtain
            // to obtain the barycentric coordinates [lambda1,lambda2] of r in the triangle.

            const double lambda1 = d11*d1+d12*d2;
            const double lambda2 = d12*d1+d22*d2;
            const double lambda0 = 1.0-lambda1-lambda2;

            return Vect3(lambda0,lambda1,lambda2);
        }

    private:

        const Triangle& triangle;

        Vect3  edge10,edge20;
        double d11,d12,d22; // Coefficients of the symmetric matrix used to find the barycentric coordinates.
    };
}
