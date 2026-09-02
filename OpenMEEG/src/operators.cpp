// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include "operators.h"
#include "monopole.h"
#include "dipole.h"

namespace OpenMEEG {

    // General routine for applying Details::operatorFerguson (see this function for further comments)
    // to an entire mesh, and storing coordinates of the output in a Matrix.

    void operatorFerguson(const Vect3& x,const Mesh& m,Matrix& mat,const unsigned& offsetI,const double coeff) {
        ThreadException e;
        #pragma omp parallel for
        #if defined NO_OPENMP || defined OPENMP_RANGEFOR
        for (const auto& vertexp : m.vertices()) {
        #elif defined OPENMP_ITERATOR
        for (auto vit=m.vertices().begin(); vit<m.vertices().end(); ++vit) {
            const Vertex* vertexp = *vit;
        #else
        for (int i=0; i<static_cast<int>(m.vertices().size()); ++i) {
            const Vertex* vertexp = *(m.vertices().begin()+i);
        #endif
            e.run([&]() {
                const unsigned vindex = vertexp->index();
                Vect3 v = Details::operatorFerguson(x,*vertexp,m);
                mat(offsetI+0,vindex) += v.x()*coeff;
                mat(offsetI+1,vindex) += v.y()*coeff;
                mat(offsetI+2,vindex) += v.z()*coeff;
            });
        }
        e.rethrow();
    }
}
