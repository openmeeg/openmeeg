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

/*! \file
    \brief file containing the integral operators
*/
#pragma once

#include <iostream>

#include <vector.h>
#include <matrix.h>
#include <symmatrix.h>
#include <sparse_matrix.h>
#include <geometry.h>
#include <integrator.h>
#include <analytics.h>

namespace OpenMEEG {

    // #define ADAPT_LHS

    // T can be a Matrix or SymMatrix
    void operatorSinternal(const Mesh& , Matrix& , const Vertices&, const double& );
    void operatorDinternal(const Mesh& , Matrix& , const Vertices&, const double& );
    void operatorFerguson(const Vect3& , const Mesh& , Matrix& , const unsigned&, const double&);
    void operatorDipolePotDer(const Vect3& , const Vect3& , const Mesh& , Vector&, const double&, const unsigned, const bool);
    void operatorDipolePot   (const Vect3& , const Vect3& , const Mesh& , Vector&, const double&, const unsigned, const bool);

    template <typename T>
    inline void _operatorD(const Triangle& T1,const Triangle& T2,T& mat,const double& coeff,const unsigned gauss_order) {
        //this version of _operatorD add in the Matrix the contribution of T2 on T1
        // for all the P1 functions it gets involved
        // consider varying order of quadrature with the distance between T1 and T2
        analyticD3 analyD(T2);

    #ifdef ADAPT_LHS
        AdaptiveIntegrator<Vect3, analyticD3> gauss(0.005);
        gauss.setOrder(gauss_order);
        Vect3 total = gauss.integrate(analyD, T1);
    #else
        STATIC_OMP Integrator<Vect3, analyticD3> gauss(gauss_order);
        Vect3 total = gauss.integrate(analyD, T1);
    #endif //ADAPT_LHS

        for (unsigned i = 0; i < 3; ++i)
            mat(T1.index(), T2(i).index()) += total(i) * coeff;
    }

    inline void _operatorDinternal(const Triangle& T2,const Vertex& P,Matrix & mat,const double& coeff) {
        analyticD3 analyD(T2);

        Vect3 total = analyD.f(P);

        for (unsigned i=0;i<3;++i)
            mat(P.index(), T2(i).index()) += total(i) * coeff;
    }

    inline double _operatorS(const Triangle& T1,const Triangle& T2,const unsigned gauss_order) {
        STATIC_OMP Triangle *oldT = 0;
        STATIC_OMP analyticS analyS;

        if ( oldT != &T1 ) { // a few computations are needed only when changing triangle T1
            oldT = (Triangle*)&T1;
            analyS.init(T1);
        }
    #ifdef ADAPT_LHS
        AdaptiveIntegrator<double, analyticS> gauss(0.005);
        gauss.setOrder(gauss_order);
        return gauss.integrate(analyS, T2);
    #else
        STATIC_OMP Integrator<double, analyticS> gauss;
        gauss.setOrder(gauss_order);
        return gauss.integrate(analyS, T2);
    #endif //ADAPT_LHS
    }

    inline double _operatorSinternal(const Triangle& T,const Vertex& P) {
        static analyticS analyS;
        analyS.init(T);
        return analyS.f(P);
    }

    template <typename T>
    inline double _operatorN(const Vertex& V1,const Vertex& V2,const Mesh& m1,const Mesh& m2,const T& mat) {

        const Mesh::VectPTriangle& trgs1 = m1.get_triangles_for_vertex(V1);
        const Mesh::VectPTriangle& trgs2 = m2.get_triangles_for_vertex(V2);

        double result = 0.0;
        for (Mesh::VectPTriangle::const_iterator tit1 = trgs1.begin(); tit1 != trgs1.end(); ++tit1 ) {
            for (Mesh::VectPTriangle::const_iterator tit2 = trgs2.begin(); tit2 != trgs2.end(); ++tit2 ) {

                // In the second case, we here divided (precalculated) operatorS by the product of areas.

                const double Iqr = (m1.current_barrier() || m2.current_barrier()) ?
                     mat((*tit1)->index()-m1.begin()->index(),(*tit2)->index()-m2.begin()->index()) :
                     mat((*tit1)->index(),(*tit2)->index())/((*tit1)->area()*(*tit2)->area());

                Vect3 CB1 = (*tit1)->next(V1)-(*tit1)->prev(V1);
                Vect3 CB2 = (*tit2)->next(V2)-(*tit2)->prev(V2);

                const double value = (CB1*CB2)*Iqr;

                // if it is the same shared vertex

                result -=  (((&m1!=&m2) && (V1==V2)) ? 0.5 : 0.25)*value;
            }
        }
        return result;
    }

    inline double _operatorP1P0(const Triangle& T2,const Vertex& V1) {
        double result = 0.;
        if  (T2.contains(V1))
            result = T2.area()/3.0;
        return result;
    }

    template <typename T>
    void operatorN(const Mesh& m1,const Mesh& m2,T& mat,const double& coeff,const unsigned gauss_order) {
        // This function has the following arguments:
        //    the 2 interacting meshes
        //    the storage Matrix for the result
        //    the coefficient to be appleid to each matrix element (depending on conductivities, ...)
        //    the gauss order parameter (for adaptive integration)

        std::cout << "OPERATOR N ... (arg : mesh " << m1.name() << " , mesh " << m2.name() << " )" << std::endl;

        unsigned i = 0; // for the PROGRESSBAR
        if (&m1==&m2) {
            if (m1.current_barrier()) {
                // we thus precompute operator S divided by the product of triangles area.

                SymMatrix matS(m1.nb_triangles());
                for (Mesh::const_iterator tit1=m1.begin();tit1!=m1.end();++tit1) {
                    PROGRESSBAR(i++, m1.nb_triangles());
                    #pragma omp parallel for
                    #ifndef OPENMP_3_0
                    for (int i2=tit1-m1.begin();i2<m1.size();++i2) {
                        const Mesh::const_iterator tit2 = m1.begin()+i2;
                    #else
                    for (Mesh::const_iterator tit2=tit1;tit2<m1.end();++tit2) {
                    #endif
                        matS(tit1->index() - m1.begin()->index(), tit2->index() - m1.begin()->index()) = _operatorS(*tit1, *tit2, gauss_order) / ( tit1->area() * tit2->area());
                    }
                }
                i = 0 ;
                for (Mesh::const_vertex_iterator vit1 = m1.vertex_begin(); vit1 != m1.vertex_end(); ++vit1) {
                    PROGRESSBAR(i++, m1.nb_vertices());
                    #pragma omp parallel for
                    #ifndef OPENMP_3_0
                    for (int i2=vit1-m1.vertex_begin();i2<m1.vertex_size();++i2) {
                        const Mesh::const_vertex_iterator vit2 = m1.vertex_begin()+i2;
                    #else
                    for (Mesh::const_vertex_iterator vit2=vit1;vit2<m1.vertex_end();++vit2) {
                    #endif
                        mat((*vit1)->index(), (*vit2)->index()) += _operatorN(**vit1, **vit2, m1, m1, matS) * coeff;
                    }
                }
            } else {
                for (Mesh::const_vertex_iterator vit1 = m1.vertex_begin(); vit1 != m1.vertex_end(); ++vit1) {
                    PROGRESSBAR(i++, m1.nb_vertices());
                    #pragma omp parallel for
                    #ifndef OPENMP_3_0
                    for (int i2=0;i2<=vit1-m1.vertex_begin();++i2) {
                        const Mesh::const_vertex_iterator vit2 = m1.vertex_begin()+i2;
                    #else
                    for (Mesh::const_vertex_iterator vit2=m1.vertex_begin();vit2<=vit1;++vit2) {
                    #endif
                        mat((*vit1)->index(), (*vit2)->index()) += _operatorN(**vit1, **vit2, m1, m1, mat) * coeff;
                    }
                }
            }
        } else {
            if ( m1.current_barrier() || m2.current_barrier() ) {
                // we thus precompute operator S divided by the product of triangles area.
                Matrix matS(m1.nb_triangles(), m2.nb_triangles());
                for (Mesh::const_iterator tit1 = m1.begin(); tit1 != m1.end(); ++tit1) {
                    PROGRESSBAR(i++, m1.nb_triangles());
                    #pragma omp parallel for
                    #ifndef OPENMP_3_0
                    for (int i2=0;i2<m2.size();++i2) {
                        const Mesh::const_iterator tit2 = m2.begin()+i2;
                    #else
                    for (Mesh::const_iterator tit2=m2.begin();tit2<m2.end();++tit2) {
                    #endif
                        matS(tit1->index()-m1.begin()->index(),tit2->index()-m2.begin()->index()) = _operatorS(*tit1,*tit2,gauss_order)/(tit1->area()*tit2->area());
                    }
                }
                i = 0 ;
                for (Mesh::const_vertex_iterator vit1 = m1.vertex_begin();vit1!=m1.vertex_end();++vit1) {
                    PROGRESSBAR(i++,m1.nb_vertices());
                    #pragma omp parallel for
                    #ifndef OPENMP_3_0
                    for (int i2=0;i2<m2.vertex_size();++i2) {
                        const Mesh::const_vertex_iterator vit2 = m2.vertex_begin()+i2;
                    #else
                    for (Mesh::const_vertex_iterator vit2=m2.vertex_begin();vit2<m2.vertex_end();++vit2) {
                    #endif
                        mat((*vit1)->index(),(*vit2)->index()) += _operatorN(**vit1,**vit2,m1,m2,matS)*coeff;
                    }
                }
            } else {
                //  This loop is exactly the same as the one just above with just matS -> mat (factorize).
                for (Mesh::const_vertex_iterator vit1 = m1.vertex_begin();vit1!=m1.vertex_end();++vit1) {
                    PROGRESSBAR(i++,m1.nb_vertices());
                    #pragma omp parallel for
                    #ifndef OPENMP_3_0
                    for (int i2=0;i2<m2.vertex_size();++i2) {
                        const Mesh::const_vertex_iterator vit2 = m2.vertex_begin()+i2;
                    #else
                    for (Mesh::const_vertex_iterator vit2=m2.vertex_begin();vit2<m2.vertex_end();++vit2) {
                    #endif
                        mat((*vit1)->index(),(*vit2)->index()) += _operatorN(**vit1,**vit2,m1,m2,mat)*coeff;
                    }
                }
            }
        }
    }

    template <typename T>
    void operatorS(const Mesh& m1, const Mesh& m2, T& mat, const double& coeff, const unsigned gauss_order)
    {
        // This function has the following arguments:
        //    the 2 interacting meshes
        //    the storage Matrix for the result
        //    the coefficient to be appleid to each matrix element (depending on conductivities, ...)
        //    the gauss order parameter (for adaptive integration)

        std::cout << "OPERATOR S ... (arg : mesh " << m1.name() << " , mesh " << m2.name() << " )" << std::endl;

        unsigned i = 0; // for the PROGRESSBAR
        // The operator S is given by Sij=\Int G*PSI(I, i)*Psi(J, j) with
        // PSI(A, a) is a P0 test function on layer A and triangle a
        if (&m1==&m2) {
            for (Mesh::const_iterator tit1=m1.begin();tit1!=m1.end();++tit1) {
                PROGRESSBAR(i++,m1.nb_triangles());
                #pragma omp parallel for
                #ifndef OPENMP_3_0
                for (int i2=tit1-m1.begin();i2<m1.size();++i2) {
                    const Mesh::const_iterator tit2 = m1.begin()+i2;
                #else
                for (Mesh::const_iterator tit2=tit1;tit2<m1.end();++tit2) {
                #endif
                    mat(tit1->index(),tit2->index()) = _operatorS(*tit1,*tit2,gauss_order)*coeff;
                }
            }
        } else {
            // TODO check the symmetry of _operatorS. 
            // if we invert tit1 with tit2: results in HeadMat differs at 4.e-5 which is too big.
            // using ADAPT_LHS with tolerance at 0.000005 (for _opS) drops this at 6.e-6. (but increase the computation time)
            for (Mesh::const_iterator tit1=m1.begin();tit1!=m1.end();++tit1) {
                PROGRESSBAR(i++,m1.nb_triangles());
                #pragma omp parallel for
                #ifndef OPENMP_3_0
                for (int i2=0;i2<m2.size();++i2) {
                    const Mesh::const_iterator tit2 = m2.begin()+i2;
                #else
                for (Mesh::const_iterator tit2=m2.begin();tit2<m2.end();++tit2) {
                #endif
                    mat(tit1->index(),tit2->index()) = _operatorS(*tit1,*tit2,gauss_order)*coeff;
                }
            }
        }
    }

    template <typename T>
    void operatorD(const Mesh& m1, const Mesh& m2, T& mat, const double& coeff, const unsigned gauss_order) {
        // This function (OPTIMIZED VERSION) has the following arguments:
        //    the 2 interacting meshes
        //    the storage Matrix for the result
        //    the coefficient to be appleid to each matrix element (depending on conductivities, ...)
        //    the gauss order parameter (for adaptive integration)

        //In this version of the function, in order to skip multiple computations of the same quantities
        //    loops are run over the triangles but the Matrix cannot be filled in this function anymore
        //    That's why the filling is done is function _operatorD
        //

        unsigned i = 0; // for the PROGRESSBAR
        #pragma omp parallel for
        #ifndef OPENMP_3_0
        for (int i1=0; i1 < m1.size(); ++i1) {
            const Mesh::const_iterator tit1 = m1.begin()+i1;
        #else
        for (Mesh::const_iterator tit1 = m1.begin(); tit1 < m1.end(); ++tit1) {
        #endif
            PROGRESSBAR(i++, m1.nb_triangles());
            for (Mesh::const_iterator tit2=m2.begin();tit2!=m2.end();++tit2)
                _operatorD(*tit1,*tit2,mat,coeff,gauss_order);
        }
    }

    template <typename T>
    void operatorD(const Mesh& m1,const Mesh& m2, T& mat,const double& coeff,const unsigned gauss_order,const bool star) {
        // This function (OPTIMIZED VERSION) has the following arguments:
        //    the 2 interacting meshes
        //    the storage Matrix for the result
        //    the coefficient to be appleid to each matrix element (depending on conductivities, ...)
        //    the gauss order parameter (for adaptive integration)
        //    an optional star parameter, which denotes the adjoint of the operator

        std::cout << "OPERATOR D" << ((star) ? "*" : " ") << "... (arg : mesh " << m1.name() << " , mesh " << m2.name() << " )" << std::endl;

        if (star) {
            operatorD(m2, m1, mat, coeff, gauss_order);
        } else {
            operatorD(m1, m2, mat, coeff, gauss_order);
        }
    }

    template <typename T>
    void operatorP1P0(const Mesh& m, T& mat,const double& coeff) {
        // This time mat(i, j)+= ... the Matrix is incremented by the P1P0 operator
        std::cout << "OPERATOR P1P0... (arg : mesh " << m.name() << " )" << std::endl;
        for (Mesh::const_iterator tit = m.begin(); tit != m.end(); ++tit)
            for (Triangle::const_iterator pit = tit->begin(); pit != tit->end(); ++pit)
                mat(tit->index(), (*pit)->index()) += _operatorP1P0(*tit, **pit) * coeff;
    }

    inline Vect3 _operatorFerguson(const Vect3& x,const Vertex& V1,const Mesh& m) {
        STATIC_OMP Vect3 result;
        STATIC_OMP analyticS analyS;

        result.x() = 0.0;
        result.y() = 0.0;
        result.z() = 0.0;

        //loop over triangles of which V1 is a vertex
        const Mesh::VectPTriangle& trgs = m.get_triangles_for_vertex(V1);

        for (Mesh::VectPTriangle::const_iterator tit=trgs.begin();tit!=trgs.end();++tit) {

            const Triangle& T1 = **tit;

            // A1 , B1  are the two opposite vertices to V1 (triangle A1, B1, V1)
            Vect3 A1   = T1.next(V1);
            Vect3 B1   = T1.prev(V1);
            Vect3 A1B1 = (A1 - B1) * (0.5 / T1.area());
            
            analyS.init(V1, A1, B1);
            const double opS = analyS.f(x);

            result += (A1B1 * opS);
        }
        return result;
    }
}
