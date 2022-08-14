// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <triangle.h>
#include <Triangle_triangle_intersection.h>

namespace OpenMEEG {

    bool Triangle::intersects(const Triangle& triangle) const {
        double* p1 = const_cast<double*>(static_cast<const double*>(vertex(0)));
        double* q1 = const_cast<double*>(static_cast<const double*>(vertex(1)));
        double* r1 = const_cast<double*>(static_cast<const double*>(vertex(2)));
        double* p2 = const_cast<double*>(static_cast<const double*>(triangle.vertex(0)));
        double* q2 = const_cast<double*>(static_cast<const double*>(triangle.vertex(1)));
        double* r2 = const_cast<double*>(static_cast<const double*>(triangle.vertex(2)));
        return tri_tri_overlap_test_3d(p1,q1,r1,p2,q2,r2);
    }
}
