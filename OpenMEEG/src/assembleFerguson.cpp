// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <operators.h>
#include <progressbar.h>
#include <constants.h>

namespace OpenMEEG {

    // geom = geometry
    // mat  = storage for Ferguson Matrix
    // pts  = where the magnetic field is to be computed

    void assemble_ferguson(const Geometry& geom,Matrix& mat,const Matrix& pts) {

        // Computation of blocks of Ferguson's Matrix

        const unsigned n = pts.nlin();
        ProgressBar pb(geom.meshes().size()*n);
        for (const auto& mesh : geom.meshes()) {
            const double coeff = MagFactor*geom.conductivity_jump(mesh);
            for (unsigned i=0,index=0; i<n; ++i,index+=3,++pb) {
                const Vect3 p(pts(i,0),pts(i,1),pts(i,2));
                operatorFerguson(p,mesh,mat,index,coeff);
            }
        }
    }
}
