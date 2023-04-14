// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <matrix.h>
#include <vector.h>
#include <vect3.h>

namespace OpenMEEG {

    class OPENMEEG_EXPORT Monopole {
    public:

        Monopole(const Vector& V): r0(V(0),V(1),V(2)),q(V(3))   { }
        Monopole(const Vect3& pos,const double charge): r0(pos),q(charge) { }

        Monopole(const unsigned i,const Matrix& M): r0(M(i,0),M(i,1),M(i,2)),q(M(i,3)) { }

        double potential(const Vect3& r) const {

            // V = q/||r-r0||

            const Vect3& x    = r-r0;
            const double nrm2 = x.norm2();
            return q/nrm2;
        }

        const Vect3& position() const { return r0; }
        double charge()         const { return q;  }

    private:

        const Vect3 r0;
        const double q;
    };
}
