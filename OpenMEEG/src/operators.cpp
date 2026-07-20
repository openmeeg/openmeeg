// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <operators.h>

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
            e.Run([&](){
                const unsigned vindex = vertexp->index();
                Vect3 v = Details::operatorFerguson(x,*vertexp,m);
                mat(offsetI+0,vindex) += v.x()*coeff;
                mat(offsetI+1,vindex) += v.y()*coeff;
                mat(offsetI+2,vindex) += v.z()*coeff;
            });
        }
        e.Rethrow();
    }

    void operatorMonopolePotDer(const Monopole& monopole,const Mesh& m,Vector& rhs,const double coeff,const Integrator& integrator) {
        ThreadException e;
        #pragma omp parallel for
        #if defined NO_OPENMP || defined OPENMP_RANGEFOR
        for (const auto& triangle : m.triangles()) {
        #elif defined OPENMP_ITERATOR
        for (Triangles::const_iterator tit=m.triangles().begin(); tit<m.triangles().end(); ++tit) {
            const Triangle& triangle = *tit;
        #else
        for (int i=0; i<static_cast<int>(m.triangles().size()); ++i) {
            const Triangle& triangle = *(m.triangles().begin()+i);
        #endif
            e.Run([&](){
                const analyticMonopolePotDer anaDPD(monopole,triangle);
                const auto monopder = [&](const Vect3& r) { return anaDPD.f(r); };

                const Vect3& v = integrator.integrate(monopder,triangle);
                // On clang/macOS we hit https://stackoverflow.com/questions/66362932/re-throwing-exception-from-openmp-block-with-the-main-thread-with-rcpp
                #ifndef __APPLE__
                #pragma omp critical
                #endif
                {
                    for (unsigned j=0; j<3; ++j)
                        rhs(triangle.vertex(j).index()) += v(j)*coeff;
                }
            });
        }
        e.Rethrow();
    }

    void operatorMonopolePot(const Monopole& monopole,const Mesh& m,Vector& rhs,const double coeff,const Integrator& integrator) {
        const auto& monoppot = [&monopole](const Vect3& r) { return monopole.potential(r); };
        ThreadException e;
        #pragma omp parallel for
        #if defined NO_OPENMP || defined OPENMP_RANGEFOR
        for (const auto& triangle : m.triangles()) {
        #elif defined OPENMP_ITERATOR
        for (Triangles::const_iterator tit=m.triangles().begin(); tit<m.triangles().end(); ++tit) {
            const Triangle& triangle = *tit;
        #else
        for (int i=0; i<static_cast<int>(m.triangles().size()); ++i) {
            const Triangle& triangle = *(m.triangles().begin()+i);
        #endif
            e.Run([&](){
                const double d = integrator.integrate(monoppot,triangle);
                rhs(triangle.index()) += d*coeff;
            });
        }
        e.Rethrow();
    }

    void operatorDipolePotDer(const Dipole& dipole,const Mesh& m,Vector& rhs,const double coeff,const Integrator& integrator) {
        ThreadException e;
        #pragma omp parallel for
        #if defined NO_OPENMP || defined OPENMP_RANGEFOR
        for (const auto& triangle : m.triangles()) {
        #elif defined OPENMP_ITERATOR
        for (Triangles::const_iterator tit=m.triangles().begin(); tit<m.triangles().end(); ++tit) {
            const Triangle& triangle = *tit;
        #else
        for (int i=0; i<static_cast<int>(m.triangles().size()); ++i) {
            const Triangle& triangle = *(m.triangles().begin()+i);
        #endif
            e.Run([&](){
                const analyticDipPotDer anaDPD(dipole,triangle);
                const auto dipder = [&](const Vect3& r) { return anaDPD.f(r); };

                const Vect3& v = integrator.integrate(dipder,triangle);
                // On clang/macOS we hit https://stackoverflow.com/questions/66362932/re-throwing-exception-from-openmp-block-with-the-main-thread-with-rcpp
                #ifndef __APPLE__
                #pragma omp critical
                #endif
                {
                    for (unsigned j=0; j<3; ++j)
                        rhs(triangle.vertex(j).index()) += v(j)*coeff;
                }
            });
        }
        e.Rethrow();
    }

    void operatorDipolePot(const Dipole& dipole,const Mesh& m,Vector& rhs,const double coeff,const Integrator& integrator) {
        const auto& dippot = [&dipole](const Vect3& r) { return dipole.potential(r); };
        ThreadException e;
        #pragma omp parallel for
        #if defined NO_OPENMP || defined OPENMP_RANGEFOR
        for (const auto& triangle : m.triangles()) {
        #elif defined OPENMP_ITERATOR
        for (Triangles::const_iterator tit=m.triangles().begin(); tit<m.triangles().end(); ++tit) {
            const Triangle& triangle = *tit;
        #else
        for (int i=0; i<static_cast<int>(m.triangles().size()); ++i) {
            const Triangle& triangle = *(m.triangles().begin()+i);
        #endif
            e.Run([&](){
                const double d = integrator.integrate(dippot,triangle);
                rhs(triangle.index()) += d*coeff;
            });
        }
        e.Rethrow();
    }
}
