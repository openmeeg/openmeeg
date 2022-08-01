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

#include <matrix.h>
#include <vector.h>
#include <vect3.h>

namespace OpenMEEG {

    class OPENMEEG_EXPORT Dipole {
    public:

        Dipole(const Vector& V): r0(V(0),V(1),V(2)),q(V(3),V(4),V(5))   { }
        Dipole(const Vect3& pos,const Vect3& moment): r0(pos),q(moment) { }

        Dipole(const unsigned i,const Matrix& M): r0(M(i,0),M(i,1),M(i,2)),q(M(i,3),M(i,4),M(i,5)) { }

        double potential(const Vect3& r) const {

            // V = q.(r-r0)/||r-r0||^3

            const Vect3& x    = r-r0;
            const double nrm2 = x.norm2();
            return dotprod(q,x)/(nrm2*sqrt(nrm2));
        }

        const Vect3& position() const { return r0; }
        const Vect3& moment()   const { return q;  }

    private:

        const Vect3 r0;
        const Vect3 q;
    };
}
