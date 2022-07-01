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

/// \file
/// \brief File containing the integral operators.

#pragma once

#include <iostream>

#include <vector.h>
#include <matrix.h>
#include <symmatrix.h>
#include <sparse_matrix.h>
#include <geometry.h>
#include <integrator.h>
#include <analytics.h>

#include <progressbar.h>

namespace OpenMEEG {

    // TODO: Use overloading and remove the internal suffix.

    void operatorSinternal(const Mesh&,Matrix&,const Vertices&,const double&);
    void operatorDinternal(const Mesh&,Matrix&,const Vertices&,const double&);
    void operatorFerguson(const Vect3&,const Mesh&,Matrix&,const unsigned&,const double&);
    void operatorDipolePotDer(const Vect3&,const Vect3&,const Mesh&,Vector&,const double&,const unsigned,const bool);
    void operatorDipolePot(const Vect3&,const Vect3&,const Mesh&,Vector&,const double&,const unsigned,const bool);

    namespace Details {
        // #define ADAPT_LHS

        template <template <typename,typename> class Integrator>
        void operatorDipolePot(const Vect3& r0,const Vect3& q,const Mesh& m,Vector& rhs,const double& coeff,const unsigned gauss_order) {
            static analyticDipPot anaDP;

            anaDP.init(q,r0);
            Integrator<double,analyticDipPot> gauss(0.001);
            gauss->setOrder(gauss_order);

            #pragma omp parallel for
            #if defined NO_OPENMP || defined OPENMP_RANGEFOR
            for (const auto& triangle : m.triangles()) {
            #elif defined OPENMP_ITERATOR
            for (Triangles::const_iterator tit=m.triangles().begin();tit<m.triangles().end();++tit) {
                const Triangle& triangle = *tit;
            #else
            for (int i=0;i<m.triangles().size();++i) {
                const Triangle& triangle = *(m.triangles().begin()+i);
            #endif
                const double d = gauss->integrate(anaDP,triangle);
                #pragma omp critical
                rhs(triangle.index()) += d*coeff;
            }
        }

        // T can be a Matrix or SymMatrix

        template <typename T>
        inline void operatorD(const Triangle& T1,const Triangle& T2,T& mat,const double& coeff,const unsigned gauss_order) {
            //this version of operatorD add in the Matrix the contribution of T2 on T1
            // for all the P1 functions it gets involved
            // consider varying order of quadrature with the distance between T1 and T2
            analyticD3 analyD(T2);

        #ifdef ADAPT_LHS
            AdaptiveIntegrator<Vect3, analyticD3> gauss(0.005);
            gauss.setOrder(gauss_order);
        #else
            STATIC_OMP Integrator<Vect3, analyticD3> gauss(gauss_order);
        #endif
            const Vect3 total = gauss.integrate(analyD,T1);

            for (unsigned i=0; i<3; ++i)
                mat(T1.index(),T2.vertex(i).index()) += total(i)*coeff;
        }

        inline double operatorS(const analyticS& analyS,const Triangle& T2,const unsigned gauss_order) {
        #ifdef ADAPT_LHS
            AdaptiveIntegrator<double,analyticS> gauss(0.005);
        #else
            STATIC_OMP Integrator<double,analyticS> gauss;
        #endif
            gauss.setOrder(gauss_order);
            return gauss.integrate(analyS,T2);
        }

        template <typename T>
        inline double operatorN(const Vertex& V1,const Vertex& V2,const Mesh& m1,const Mesh& m2,const T& mat) {

            const bool same_shared_vertex = ((&m1!=&m2) && (V1==V2));
            const double factor = (same_shared_vertex) ? 0.5 : 0.25;

            double result = 0.0;
            for (const auto& tp1 : m1.triangles(V1)) {
                const Edge& edge1 = tp1->edge(V1);
                const Vect3& CB1 = edge1.vertex(0)-edge1.vertex(1);
                const unsigned ind1 = tp1->index()-m1.triangles().front().index();
                for (const auto& tp2 : m2.triangles(V2)) {

                    const unsigned ind2 = tp2->index()-m2.triangles().front().index();

                    // In the second case, we here divided (precalculated) operatorS by the product of areas.

                    const double Iqr = (m1.current_barrier() || m2.current_barrier()) ? mat(ind1,ind2) : mat(tp1->index(),tp2->index())/(tp1->area()*tp2->area());

                    const Edge& edge2 = tp2->edge(V2);
                    const Vect3& CB2 = edge2.vertex(0)-edge2.vertex(1);

                    result -= factor*Iqr*dotprod(CB1,CB2);
                }
            }
            return result;
        }

        inline Vect3 operatorFerguson(const Vect3& x,const Vertex& V,const Mesh& m) {
            Vect3 result;
            result = 0.0;

            //  Loop over triangles of which V is a vertex

            for (const auto& tp : m.triangles(V)) {
                const Triangle& T    = *tp;
                const Edge&     edge = T.edge(V);

                // A, B are the two opposite vertices to V (triangle A, B, V)
                const Vertex& A = edge.vertex(0);
                const Vertex& B = edge.vertex(1);
                const Vect3 AB = (A-B)*(0.5/T.area());
                
                analyticS analyS(V,A,B);
                const double opS = analyS.f(x);

                result += (AB*opS);
            }

            return result;
        }

        inline double operatorP1P0(const Triangle& T2,const Vertex& V1) {
            double result = 0.;
            if  (T2.contains(V1))
                result = T2.area()/3.0;
            return result;
        }

        template <typename T>
        void operatorD(const Mesh& m1,const Mesh& m2,T& mat,const double& coeff,const unsigned gauss_order) {
            // This function (OPTIMIZED VERSION) has the following arguments:
            //    the 2 interacting meshes
            //    the storage Matrix for the result
            //    the coefficient to be appleid to each matrix element (depending on conductivities, ...)
            //    the gauss order parameter (for adaptive integration)

            // In this version of the function, in order to skip multiple computations of the same quantities
            //    loops are run over the triangles but the Matrix cannot be filled in this function anymore
            //    That's why the filling is done is function Details::operatorD
            //

            ProgressBar pb(m1.triangles().size());
            const Triangles& m1_triangles = m1.triangles();
            #pragma omp parallel for
            #if defined NO_OPENMP || defined OPENMP_RANGEFOR
            for (const auto& triangle1 : m1_triangles) {
            #elif defined OPENMP_ITERATOR
            for (Triangles::const_iterator tit1=m1_triangles.begin();tit1<m1_triangles.end();++tit1) {
                const Triangle& triangle1 = *tit1;
            #else
            for (int i1=0; i1 < m1_triangles.size(); ++i1) {
                const Triangle& triangle1 = *(m1_triangles.begin()+i1);
            #endif
                for (const auto& triangle2 : m2.triangles())
                    Details::operatorD(triangle1,triangle2,mat,coeff,gauss_order);
                ++pb;
            }
        }
    }

    template <typename T>
    void operatorN(const Mesh& m1,const Mesh& m2,T& mat,const double& coeff,const unsigned gauss_order) {
        // This function has the following arguments:
        //    the 2 interacting meshes
        //    the storage Matrix for the result
        //    the coefficient to be applied to each matrix element (depending on conductivities, ...)
        //    the gauss order parameter (for adaptive integration)

        std::cout << "OPERATOR N ... (arg : mesh " << m1.name() << " , mesh " << m2.name() << " )" << std::endl;

        if (&m1==&m2) {
            auto NUpdate = [&](const Mesh& m,const auto& M) {
                ProgressBar pb(m1.vertices().size());
                for (auto vit1=m.vertices().begin();vit1!=m.vertices().end();++vit1) {
                    #pragma omp parallel for
                    #if defined NO_OPENMP || defined OPENMP_ITERATOR
                    for (auto vit2=vit1;vit2<m.vertices().end();++vit2) {
                    #else
                    for (int i2=0;i2<=vit1-m1.vertices().begin();++i2) {
                        const auto vit2 = m1.vertices().begin()+i2;
                    #endif
                        mat((*vit1)->index(),(*vit2)->index()) += Details::operatorN(**vit1,**vit2,m,m,M)*coeff;
                    }
                    ++pb;
                }
            };

            if (m1.current_barrier()) {
                // Precompute operator S divided by the product of triangles area.

                ProgressBar pb(m1.triangles().size());
                SymMatrix matS(m1.triangles().size());
                for (Triangles::const_iterator tit1=m1.triangles().begin();tit1!=m1.triangles().end();++tit1) {
                    const analyticS analyS(*tit1);
                    const unsigned ind1 = tit1->index()-m1.triangles().front().index();
                    #pragma omp parallel for
                    #if defined NO_OPENMP || defined OPENMP_ITERATOR
                    for (Triangles::const_iterator tit2=tit1;tit2<m1.triangles().end();++tit2) {
                    #else
                    for (int i2=tit1-m1.triangles().begin();i2<m1.triangles().size();++i2) {
                        const Triangles::const_iterator tit2 = m1.triangles().begin()+i2;
                    #endif
                        const unsigned ind2 = tit2->index()-m2.triangles().front().index();
                        matS(ind1,ind2) = Details::operatorS(analyS,*tit2,gauss_order)/(tit1->area()*tit2->area());
                    }
                    ++pb;
                }
                NUpdate(m1,matS);
            } else {
                NUpdate(m1,mat);
            }
        } else {
            auto NUpdate = [&](const Mesh& m1,const Mesh& m2,const auto& M) {
                ProgressBar pb(m1.vertices().size());
                const VerticesRefs& v2 = m2.vertices();
                for (const auto& vertex1 : m1.vertices()) {
                    #pragma omp parallel for
                    #if defined NO_OPENMP || defined OPENMP_RANGEFOR
                    for (const auto& vertex2 : v2) {
                    #elif defined OPENMP_ITERATOR
                    for (auto vit2=v2.begin();vit2<v2.end();++vit2) {
                        const Vertex* vertex2 = *vit2;
                    #else
                    for (int i2=0;i2<v2.size();++i2) {
                        const Vertex* vertex2 = *(v2.begin()+i2);
                    #endif
                        mat(vertex1->index(),vertex2->index()) += Details::operatorN(*vertex1,*vertex2,m1,m2,M)*coeff;
                    }
                    ++pb;
                }
            };

            if (m1.current_barrier() || m2.current_barrier()) {
                // Precompute operator S divided by the product of triangles area.
                Matrix matS(m1.triangles().size(),m2.triangles().size());
                ProgressBar pb(m1.triangles().size());
                unsigned i = 0;
                for (const auto& triangle1 : m1.triangles()) {
                    const analyticS analyS(triangle1);
                    const unsigned ind1 = triangle1.index()-m1.triangles().front().index();
                    const Triangles& m2_triangles = m2.triangles();
                    #pragma omp parallel for
                    #if defined NO_OPENMP || defined OPENMP_RANGEFOR
                    for (const auto& triangle2 : m2_triangles) {
                    #elif defined OPENMP_ITERATOR
                    for (Triangles::const_iterator tit2=m2_triangles.begin();tit2<m2_triangles.end();++tit2) {
                        const Triangle& triangle2 = *tit2;
                    #else
                    for (int i2=0;i2<m2_triangles.size();++i2) {
                        const Triangle& triangle2 = *(m2_triangles.begin()+i2);
                    #endif
                        const unsigned ind2 = triangle2.index()-m2_triangles.front().index();
                        matS(ind1,ind2) = Details::operatorS(analyS,triangle2,gauss_order)/(triangle1.area()*triangle2.area());
                    }
                    ++pb;
                }
                NUpdate(m1,m2,matS);
            } else {
                NUpdate(m1,m2,mat);
            }
        }
    }

    template <typename T>
    void operatorS(const Mesh& m1,const Mesh& m2,T& mat,const double& coeff,const unsigned gauss_order) {

        // This function has the following arguments:
        //    the 2 interacting meshes
        //    the storage Matrix for the result
        //    the coefficient to be applied to each matrix element (depending on conductivities, ...)
        //    the gauss order parameter (for adaptive integration)

        std::cout << "OPERATOR S ... (arg : mesh " << m1.name() << " , mesh " << m2.name() << " )" << std::endl;

        // The operator S is given by Sij=\Int G*PSI(I, i)*Psi(J, j) with
        // PSI(A, a) is a P0 test function on layer A and triangle a
        if (&m1==&m2) {
            ProgressBar pb(m1.triangles().size());
            for (Triangles::const_iterator tit1=m1.triangles().begin(); tit1!=m1.triangles().end(); ++tit1,++pb) {
                const analyticS analyS(*tit1);
                #pragma omp parallel for
                #if defined OPENMP_ITERATOR
                for (Triangles::const_iterator tit2=tit1;tit2<m1.triangles().end();++tit2) {
                #else
                for (int i2=tit1-m1.triangles().begin();i2<m1.triangles().size();++i2) {
                    const Triangles::const_iterator tit2 = m1.triangles().begin()+i2;
                #endif
                    mat(tit1->index(),tit2->index()) = Details::operatorS(analyS,*tit2,gauss_order)*coeff;
                }
            }
        } else {
            // TODO check the symmetry of Details::operatorS. 
            // if we invert tit1 with tit2: results in HeadMat differs at 4.e-5 which is too big.
            // using ADAPT_LHS with tolerance at 0.000005 (for Details::opS) drops this at 6.e-6. (but increase the computation time)

            ProgressBar pb(m1.triangles().size());
            for (const auto& triangle1 : m1.triangles()) {
                const analyticS analyS(triangle1);
                const Triangles& m2_triangles = m2.triangles();
                #pragma omp parallel for
                #if defined NO_OPENMP || defined OPENMP_RANGEFOR
                for (const auto& triangle2 : m2_triangles) {
                #elif defined OPENMP_ITERATOR
                for (Triangles::const_iterator tit2=m2_triangles.begin();tit2<m2_triangles.end();++tit2) {
                    const Triangle& triangle2 = *tit2;
                #else
                for (int i2=0;i2<m2_triangles.size();++i2) {
                    const Triangle& triangle2 = *(m2_triangles.begin()+i2);
                #endif
                    mat(triangle1.index(),triangle2.index()) = Details::operatorS(analyS,triangle2,gauss_order)*coeff;
                }
                ++pb;
            }
        }
    }

    template <typename T>
    void operatorD(const Mesh& m1,const Mesh& m2,T& mat,const double& coeff,const unsigned gauss_order) {
        // This function (OPTIMIZED VERSION) has the following arguments:
        //    the 2 interacting meshes
        //    the storage Matrix for the result
        //    the coefficient to be appleid to each matrix element (depending on conductivities, ...)
        //    the gauss order parameter (for adaptive integration)

        std::cout << "OPERATOR D... (arg : mesh " << m1.name() << " , mesh " << m2.name() << " )" << std::endl;
        Details::operatorD(m1,m2,mat,coeff,gauss_order);
    }

    template <typename T>
    void operatorDstar(const Mesh& m1,const Mesh& m2,T& mat,const double& coeff,const unsigned gauss_order) {
        // This function (OPTIMIZED VERSION) has the following arguments:
        //    the 2 interacting meshes
        //    the storage Matrix for the result
        //    the coefficient to be appleid to each matrix element (depending on conductivities, ...)
        //    the gauss order parameter (for adaptive integration)

        std::cout << "OPERATOR D*... (arg : mesh " << m1.name() << " , mesh " << m2.name() << ')' << std::endl;
        Details::operatorD(m2,m1,mat,coeff,gauss_order);
    }

    template <typename T>
    void operatorP1P0(const Mesh& m,T& mat,const double& coeff) {
        // This time mat(i, j)+= ... the Matrix is incremented by the P1P0 operator
        std::cout << "OPERATOR P1P0... (arg : mesh " << m.name() << " )" << std::endl;
        for (const auto& triangle : m.triangles())
            for (const auto& vertex : triangle)
                mat(triangle.index(),vertex->index()) += Details::operatorP1P0(triangle,*vertex)*coeff;
    }
}
