// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

/// \file
/// \brief File containing the integral operators.

#include <iostream>
#include <iomanip>

#include <vector.h>
#include <matrix.h>
#include <symmatrix.h>
#include <block.h>
#include <geometry.h>
#include <integrator.h>
#include <analytics.h>
#include <constants.h>

#include <logger.h>
#include <progressbar.h>

namespace OpenMEEG {

    void operatorFerguson(const Vect3&,const Mesh&,Matrix&,const unsigned&,const double);

    template <typename Source>
    Vector operatorPotential(const unsigned size,const Source& source,const Mesh& m,const Integrator& integrator) {
        Vector res(size);
        res = 0.0;
        const auto& potential = [&source](const Point3D& r) { return source.potential(r); };
        ThreadException e;
        #pragma omp parallel for
        #if defined NO_OPENMP || defined OPENMP_RANGEFOR
        for (const auto& triangle : m.triangles()) {
        #else
        for (int i=0; i<static_cast<int>(m.triangles().size()); ++i) {
            const Triangle& triangle = *(m.triangles().begin()+i);
        #endif
            e.run([&]() {
                res(triangle.index()) += integrator.integrate(potential,triangle);
            });
        }
        e.rethrow();
        return res;
    }

    template <typename Source>
    Vector operatorPotentialDerivative(const unsigned size,const Source& source,const Mesh& m,const Integrator& integrator) {
        Vector res(size);
        res = 0.0;
        ThreadException e;
        #pragma omp parallel for
        #if defined NO_OPENMP || defined OPENMP_RANGEFOR
        for (const auto& triangle : m.triangles()) {
        #else
        for (int i=0; i<static_cast<int>(m.triangles().size()); ++i) {
            const Triangle& triangle = *(m.triangles().begin()+i);
        #endif
            e.run([&]() {
                const PotentialDerivative<Source> potential_derivative(source,triangle);
                const auto pot_derivative = [&](const Point3D& r) { return potential_derivative(r); };
                const Vect3& v = integrator.integrate(pot_derivative,triangle);

                for (unsigned j=0; j<3; ++j) {
                    double& r = res(triangle.vertex(j).index());
                    #pragma omp atomic
                    r += v(j);
                }
            });
        }
        e.rethrow();
        return res;
    }

    namespace Details {

        inline Vect3 operatorFerguson(const Vect3& x,const Vertex& V,const Mesh& m) {
            Vect3 result;

            //  Loop over triangles of which V is a vertex

            result = 0.0;
            for (const auto& tp : m.triangles(V)) {
                const Triangle& T    = *tp;
                const Edge&     edge = T.edge(V);

                // A, B are the two opposite vertices to V in T (triangle A, B, V)

                const Vertex& A = edge.vertex(0);
                const Vertex& B = edge.vertex(1);

                const OperatorS S(T);

                result += (B-A)*S(x)/(2*T.area());
            }

            return result;
        }
    }

    class BlockWriterBase {
    public:

        BlockWriterBase(const Integrator& intg): integrator(intg) { }

        void message(const char* op_name,const Mesh& mesh) const {
            log_stream(INFORMATION) << std::endl
                                    << "OPERATOR " << std::left << std::setw(2) << op_name
                                    << "... (arg : mesh " << mesh.name() << " )" << std::endl;
        }

        void message(const char* op_name,const Mesh& mesh1,const Mesh& mesh2) const {
            log_stream(INFORMATION) << "OPERATOR " << std::left << std::setw(2) << op_name
                                    << "... (arg : mesh " << mesh1.name() << " , mesh " << mesh2.name() << " )"
                                    << std::endl;
        }

    protected:

        template <typename DMatrix>
        void D(const Triangles& triangles1,const Triangles& triangles2,const double coeff,DMatrix& Dmat) const {
            // This function (OPTIMIZED VERSION) has the following arguments:
            //    - 2 interacting meshes (as sets of triangles).
            //    - coefficient to be applied to each matrix element (depending on conductivities, ...)
            //    - storage Matrix for the result.

            ThreadException e;
            ProgressBar pb(triangles1.size());
            #pragma omp parallel for
            #if defined NO_OPENMP || defined OPENMP_RANGEFOR
            for (const auto& triangle1 : triangles1) {
            #else
            for (int i1=0; i1<static_cast<int>(triangles1.size()); ++i1) {
                const Triangle& triangle1 = *(triangles1.begin()+i1);
            #endif
                e.run([&]() {
                    for (const auto& triangle2 : triangles2) {
                        const P1SingleLayerPotential P1_single_layer_potential(triangle2);
                        const auto&  P1SLP = [&P1_single_layer_potential](const Point3D& r) { return P1_single_layer_potential(r); };
                        const Vect3& total = integrator.integrate(P1SLP,triangle1);

                        for (unsigned i=0; i<3; ++i)
                            Dmat(triangle1.index(),triangle2.vertex(i).index()) += total(i)*coeff;
                    }
                });
                ++pb;
            }
            e.rethrow();
        }

        template <typename SMatrix>
        static double N(const Vertex& V1,const Vertex& V2,const Mesh& m1,const Mesh& m2,const SMatrix& Smat) {

            double result = 0.0;
            for (const auto& tp1 : m1.triangles(V1)) {
                const Edge& edge1 = tp1->edge(V1);
                const Vect3& CB1 = edge1.vertex(0)-edge1.vertex(1);
                for (const auto& tp2 : m2.triangles(V2)) {
                    const Edge& edge2 = tp2->edge(V2);
                    const Vect3& CB2 = edge2.vertex(0)-edge2.vertex(1);
                    result -= dotprod(CB1,CB2)*Smat(tp1->index(),tp2->index())/(tp1->area()*tp2->area());
                }
            }
            return result;
        }

        const Integrator integrator;
    };

    class DiagonalBlockWriter: public BlockWriterBase {

        using base = BlockWriterBase;

    public:

        DiagonalBlockWriter(const Mesh& m,const Integrator& intg): base(intg),mesh(m) { }

        bool has_S_block() const { return !mesh.is_current_barrier(); }
        bool has_D_block() const { return !mesh.is_current_barrier(); }

        template <typename SMatrix>
        void set_S_block(const double coeff,SMatrix& Smat) {
            if (has_S_block()) {
                S(coeff,Smat);
                Scoeff = coeff;
            }
        }

        template <typename SMatrix,typename NMatrix>
        void set_N_block(const double coeff,const SMatrix& Smat,NMatrix& Nmat) const {
            if (S_block_is_computed()) {
                N(coeff/Scoeff,Smat,Nmat);
            } else {
                set_N_block(coeff,Nmat);
            }
        }

        template <typename NMatrix>
        void set_N_block(const double coeff,NMatrix& Nmat) const {

            // In this case, Smat is invalid in the FullMatrix representation and not filled
            // in the BlockMatrix representation. We need to compute this matrix.

            SymMatrix Smat(mesh.triangles().size());
            auto ind = [&](const unsigned i) { return i-mesh.triangles().front().index(); };
            Block<SymMatrix> Sblock(Smat,ind,ind);
            S(1.0,Sblock);
            N(coeff,Sblock,Nmat);
        }

        template <typename DMatrix>
        void set_D_block(const double coeff,DMatrix& Dmat) const {
            if (has_D_block()) {
                base::message("D",mesh,mesh);
                base::D(mesh.triangles(),mesh.triangles(),coeff,Dmat);
            }
        }

        // There is no need to compute a Dstar block in this case as it is the transpose(D).

        template <typename T>
        void addIdentity(const double coeff,T& matrix) const {
            // The Matrix is incremented by the identity P1P0 operator
            base::message("Id",mesh);
            const double factor = coeff/3.0;
            for (const auto& triangle : mesh.triangles())
                for (const auto& vertex : triangle)
                    matrix(triangle.index(),vertex->index()) += (triangle.contains(*vertex)) ? factor*triangle.area() : 0.0;
        }

        template <typename DMatrix>
        void set_D_block(DMatrix& Dmat) const {
            base::message("D",mesh,mesh);
            base::D(mesh.triangles(),mesh.triangles(),K,Dmat);
        }

    private:

        bool S_block_is_computed() const { return Scoeff!=0.0; }

        template <typename SMatrix>
        void S(const double coeff,SMatrix& Smat) const {
            base::message("S",mesh,mesh);
            ProgressBar pb(mesh.triangles().size());

            // Operator S is given by Sij=\Int G*PSI(I,i)*Psi(J,j) with PSI(l,t) a P0 test function on layer l and triangle t.
            // When meshes are equal, optimized computation for a symmetric matrix.

            ThreadException e;
            const Triangles& triangles = mesh.triangles();
            #pragma omp parallel for
            #if defined NO_OPENMP || defined OPENMP_ITERATOR
            for (Triangles::const_iterator tit1=triangles.begin(); tit1!=triangles.end(); ++tit1) {
                const Triangle& triangle1 = *tit1;
            #else
            for (int i1=0; i1<triangles.size(); ++i1,++pb) {
                const Triangle& triangle1 = triangles[i1];
                const Triangles::const_iterator tit1 = triangles.cbegin()+i1;
            #endif
                e.run([&]() {
                    const OperatorS S(triangle1);
                    const auto& Sfunc = [&S](const Point3D& r) { return S(r); };

                    for (Triangles::const_iterator tit2=tit1; tit2!=triangles.end(); ++tit2) {
                        const Triangle& triangle2 = *tit2;
                        Smat(triangle1.index(),triangle2.index()) = base::integrator.integrate(Sfunc,triangle2)*coeff;
                    }
                });
                ++pb;
            }
            e.rethrow();
        }

        template <typename SMatrix,typename NMatrix>
        void N(const double coeff,const SMatrix& Smat,NMatrix& Nmat) const {

            base::message("N",mesh,mesh);

            ThreadException e;
            ProgressBar pb(mesh.vertices().size());

            // When meshes are equal, optimized computation for a symmetric matrix.

            const double factor = 0.25*coeff;
            #pragma omp parallel for
            #if defined NO_OPENMP || defined OPENMP_ITERATOR
            for (auto vit1=mesh.vertices().begin(); vit1!=mesh.vertices().end(); ++vit1) {
            #else
            for (int i1=0; i1!=mesh.vertices().size(); ++i1) {
                const auto vit1 = mesh.vertices().begin()+i1;
            #endif
                e.run([&]() {
                    for (auto vit2=vit1; vit2<mesh.vertices().end(); ++vit2)
                        Nmat((*vit1)->index(),(*vit2)->index()) += factor*base::N(**vit1,**vit2,mesh,mesh,Smat);
                });
                ++pb;
            }
            e.rethrow();
        }

        const Mesh&  mesh;
              double Scoeff = 0.0;
    };

    class PartialBlockWriter {
    public:

        PartialBlockWriter(const Mesh& m): mesh(m) { }

        void addD(const double coeff,const Vertices& points,Matrix& matrix) const {
            log_stream(INFORMATION) << "PARTAL OPERATOR D..." << std::endl;
            for (const auto& triangle : mesh.triangles()) {
                const P1SingleLayerPotential P1_single_layer_potential(triangle);
                for (const auto& vertex : points) {
                    const Vect3& integrals = P1_single_layer_potential(vertex);
                    for (unsigned i=0;i<3;++i)
                        matrix(vertex.index(),triangle.vertex(i).index()) += integrals(i)*coeff;
                }
            }
        }

        void S(const double coeff,const Vertices& points,Matrix& matrix) const {
            log_stream(INFORMATION) << "PARTIAL OPERATOR S..." << std::endl;
            for (const auto& triangle : mesh.triangles()) {
                const OperatorS S(triangle);
                for (const auto& vertex : points)
                    matrix(vertex.index(),triangle.index()) = coeff*S(vertex);
            }
        }

    private:

        const Mesh& mesh;
    };

    class NonDiagonalBlockWriter: public BlockWriterBase  {

        using base = BlockWriterBase;

    public:

        // This constructor takes the following arguments:
        //  - The 2 adjacent (interacting) meshes.
        //  - The integrator used to compute the operators.

        NonDiagonalBlockWriter(const Mesh& m1,const Mesh& m2,const Integrator& intg): base(intg),mesh1(m1),mesh2(m2) { }

        template <typename SMatrix>
        void set_S_block(const double coeff,SMatrix& Smat) {
            if (has_S_block()) {
                S(coeff,Smat);
                Scoeff = coeff;
            }
        }

        template <typename SMatrix,typename NMatrix>
        void set_N_block(const double coeff,const SMatrix& Smat,NMatrix& Nmat) const {
            if (S_block_is_computed()) {
                N(coeff/Scoeff,Smat,Nmat);
            } else {
                set_N_block(coeff,Nmat);
            }
        }

        template <typename NMatrix>
        void set_N_block(const double coeff,NMatrix& Nmat) const {

            // In this case, Smat is invalid in the FullMatrix representation and not filled
            // in the BlockMatrix representation. We need to compute this matrix.

            Matrix Smat(mesh1.triangles().size(),mesh2.triangles().size());
            auto ind1 = [&](const unsigned i) { return i-mesh1.triangles().front().index(); };
            auto ind2 = [&](const unsigned i) { return i-mesh2.triangles().front().index(); };
            Block<Matrix> Sblock(Smat,ind1,ind2);
            S(1.0,Sblock);
            N(coeff,Sblock,Nmat);
        }

        template <typename DMatrix>
        void set_D_block(const double coeff,DMatrix& Dmat) const {
            if (has_D_block())
                D(coeff,Dmat);
        }

        template <typename DstarMatrix>
        void set_Dstar_block(const double coeff,DstarMatrix& Dstarmat) const {
            if (has_Dstar_block())
                Dstar(coeff,Dstarmat);
        }

        // Various operators take the following arguments:
        //  - The coefficient to be applied to each matrix element (depending on conductivities, ...)
        //  - The storage Matrix for the result

        template <typename SMatrix>
        void S(const double coeff,SMatrix& Smat) const {
            base::message("S",mesh1,mesh2);
            ThreadException e;
            ProgressBar pb(mesh1.triangles().size());

            // Operator S is given by Sij=\Int G*PSI(I,i)*Psi(J,j) with PSI(l,t) a P0 test function on layer l and triangle t.

            // TODO check the symmetry of S.
            // if we invert tit1 with tit2: results in HeadMat differs at 4.e-5 which is too big.
            // using ADAPT_LHS with tolerance at 0.000005 (for S) drops this at 6.e-6 (but increase the computation time).

            #pragma omp parallel for
            #if defined NO_OPENMP || defined OPENMP_RANGEFOR
            for (const Triangle& triangle1 : mesh1.triangles()) {
            #else
            for (int i=0; i<mesh1.triangles().size(); ++i) {
                const Triangle& triangle1 = mesh1.triangles()[i];
            #endif
                e.run([&]() {
                    const OperatorS S(triangle1);
                    const auto& Sfunc = [&S](const Point3D& r) { return S(r); };

                    for (const auto& triangle2 : mesh2.triangles())
                        Smat(triangle1.index(),triangle2.index()) = base::integrator.integrate(Sfunc,triangle2)*coeff;
                });
                ++pb;
            }
            e.rethrow();
        }

        template <typename SMatrix,typename NMatrix>
        void N(const double coeff,const SMatrix& Smat,NMatrix& Nmat) const {

            base::message("N",mesh1,mesh2);

            ThreadException e;
            ProgressBar pb(mesh1.vertices().size());
            #pragma omp parallel for
            #if defined NO_OPENMP || defined OPENMP_RANGEFOR
            for (const Vertex* vertex1 : mesh1.vertices()) {
            #else
            for (int i=0; i<mesh1.vertices().size(); ++i) {
                const Vertex* vertex1 = mesh1.vertices()[i];
            #endif
                e.run([&]() {
                    for (const auto& vertex2 : mesh2.vertices()) {
                        const double factor = ((vertex1==vertex2) ? 0.5 : 0.25)*coeff;
                        Nmat(vertex1->index(),vertex2->index()) += factor*base::N(*vertex1,*vertex2,mesh1,mesh2,Smat);
                    }
                });
                ++pb;
            }
            e.rethrow();
        }

        template <typename DMatrix>
        void D(const double coeff,DMatrix& Dmat) const {
            base::message("D",mesh1,mesh2);
            base::D(mesh1.triangles(),mesh2.triangles(),coeff,Dmat);
        }

        template <typename DstarMatrix>
        void Dstar(const double coeff,DstarMatrix& Dstarmat) const {
            base::message("D*",mesh1,mesh2);
            base::D(mesh2.triangles(),mesh1.triangles(),coeff,Dstarmat);
        }

    private:

        bool has_S_block()     const { return !(mesh1.is_current_barrier() || mesh2.is_current_barrier()); }
        bool has_D_block()     const { return (!mesh1.is_current_barrier()); }
        bool has_Dstar_block() const { return (!mesh2.is_current_barrier()); }

        bool S_block_is_computed() const { return Scoeff!=0.0; }

        const Mesh&  mesh1;
        const Mesh&  mesh2;
              double Scoeff = 0.0;
    };
}
