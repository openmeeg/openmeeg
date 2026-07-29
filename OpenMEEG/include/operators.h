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
#include <symm_block_matrix.h>
#include <sparse_matrix.h>
#include <geometry.h>
#include <integrator.h>
#include <analytics.h>

#include <logger.h>
#include <progressbar.h>

namespace OpenMEEG {

    void operatorFerguson(const Vect3&,const Mesh&,Matrix&,const unsigned&,const double);
    void operatorMonopolePotDer(const Monopole&,const Mesh&,Vector&,const double,const Integrator&);
    void operatorMonopolePot(const Monopole&,const Mesh&,Vector&,const double,const Integrator&);
    void operatorDipolePotDer(const Dipole&,const Mesh&,Vector&,const double,const Integrator&);
    void operatorDipolePot(const Dipole&,const Mesh&,Vector&,const double,const Integrator&);

    namespace Details {

        inline Vect3 operatorFerguson(const Vect3& x,const Vertex& V,const Mesh& m) {
            Vect3 result;

            //  Loop over triangles of which V is a vertex

            result = 0.0;
            for (const auto& tp : m.triangles(V)) {
                const Triangle& T    = *tp;
                const Edge&     edge = T.edge(V);

                // A, B are the two opposite vertices to V (triangle A, B, V)
                const Vertex& A = edge.vertex(0);
                const Vertex& B = edge.vertex(1);
                const Vect3& AB = (A-B)/(2*T.area());

                const analyticS analyS(V,A,B);

                result += AB*analyS.f(x);
            }

            return result;
        }
    }

    class BlocksBase {
    public:

        BlocksBase(const Integrator& intg): integrator(intg) { }

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

        template <typename T>
        void D(const Triangles& triangles1,const Triangles& triangles2,const double coeff,T& mat) const {
            // This function (OPTIMIZED VERSION) has the following arguments:
            //    - 2 interacting meshes (as sets of triangles).
            //    - coefficient to be applied to each matrix element (depending on conductivities, ...)
            //    - storage Matrix for the result.

            ThreadException e;
            ProgressBar pb(triangles1.size());
            #pragma omp parallel for
            #if defined NO_OPENMP || defined OPENMP_RANGEFOR
            for (const auto& triangle1 : triangles1) {
            #elif defined OPENMP_ITERATOR
            for (Triangles::const_iterator tit1=triangles1.begin(); tit1<triangles1.end(); ++tit1) {
                const Triangle& triangle1 = *tit1;
            #else
            for (int i1=0; i1<static_cast<int>(triangles1.size()); ++i1) {
                const Triangle& triangle1 = *(triangles1.begin()+i1);
            #endif
                e.Run([&](){
                    for (const auto& triangle2 : triangles2) {
                        const analyticD3 analyD(triangle2);
                        const auto&  Dfunc = [&analyD](const Vect3& r) { return analyD.f(r); };
                        const Vect3& total = integrator.integrate(Dfunc,triangle1);

                        for (unsigned i=0; i<3; ++i)
                            mat(triangle1.index(),triangle2.vertex(i).index()) += total(i)*coeff;
                    }
                    ++pb;
                });
            }
            e.Rethrow();
        }

        // Operator N for two vertices of the same mesh.

        template <typename T>
        static double N(const Vertex& V1,const Vertex& V2,const Mesh& m,const T& matrix) {
            return N(0.25,V1,V2,m,m,matrix);
        }

        // This function assumes that the two meshes are different.

        template <typename T>
        static double N(const Vertex& V1,const Vertex& V2,const Mesh& m1,const Mesh& m2,const T& matrix) {
            const double coeff = (V1==V2) ? 0.5 : 0.25;
            return N(coeff,V1,V2,m1,m2,matrix);
        }

    private:

        template <typename T>
        static double N(const double factor,const Vertex& V1,const Vertex& V2,const Mesh& m1,const Mesh& m2,const T& matrix) {

            double result = 0.0;
            for (const auto& tp1 : m1.triangles(V1)) {
                const Edge& edge1 = tp1->edge(V1);
                const Vect3& CB1 = edge1.vertex(0)-edge1.vertex(1);
                for (const auto& tp2 : m2.triangles(V2)) {
                    const Edge& edge2 = tp2->edge(V2);
                    const Vect3& CB2 = edge2.vertex(0)-edge2.vertex(1);
                    result -= factor*dotprod(CB1,CB2)*matrix(tp1->index(),tp2->index())/(tp1->area()*tp2->area());
                }
            }
            return result;
        }

    protected:

        const Integrator integrator;
    };

    class DiagonalBlock: public BlocksBase {

        typedef BlocksBase base;

        class SymBloc: public SymMatrix {

            typedef SymMatrix base;

        public:

            SymBloc(const unsigned off,const unsigned sz): base(sz),offset(off) { }

            double& operator()(const unsigned i,const unsigned j)       { return base::operator()(i-offset,j-offset); }
            double  operator()(const unsigned i,const unsigned j) const { return base::operator()(i-offset,j-offset); }

        private:

            const unsigned offset;
        };

    public:

        DiagonalBlock(const Mesh& m,const Integrator& intg): base(intg),mesh(m) { }

        template <typename T>
        void set_S_block(const double coeff,T& matrix) {
            if (!mesh.current_barrier()) {
                S(coeff,matrix);
                Scoeff = coeff;
            }
        }

        template <typename T>
        void set_N_block(const double coeff,T& matrix) const { N(coeff,matrix); }

        template <typename T>
        void set_D_block(const double coeff,T& matrix) const {
            if (!mesh.current_barrier())
                D(coeff,matrix);
        }

        template <typename T>
        void set_Dstar_block(const double /* coeff */,T& /* matrix */) const { }

        template <typename T>
        void addIdentity(const double coeff,T& matrix) const {
            // The Matrix is incremented by the identity P1P0 operator
            base::message("Id",mesh);
            for (const auto& triangle : mesh.triangles())
                for (const auto& vertex : triangle)
                    matrix(triangle.index(),vertex->index()) += Id(triangle,*vertex)*coeff;
        }

        template <typename T>
        void S(const double coeff,T& matrix) const {
            base::message("S",mesh,mesh);
            ProgressBar pb(mesh.triangles().size());

            // Operator S is given by Sij=\Int G*PSI(I,i)*Psi(J,j) with PSI(l,t) a P0 test function on layer l and triangle t.
            // When meshes are equal, optimized computation for a symmetric matrix.

            ThreadException e;
            const Triangles& triangles = mesh.triangles();
            for (Triangles::const_iterator tit1=triangles.begin(); tit1!=triangles.end(); ++tit1,++pb) {
                const Triangle& triangle1 = *tit1;
                const analyticS analyS(triangle1);
                const auto& Sfunc = [&analyS](const Vect3& r) { return analyS.f(r); };

                #pragma omp parallel for
                #if defined NO_OPENMP || defined OPENMP_ITERATOR
                for (Triangles::const_iterator tit2=tit1; tit2!=triangles.end(); ++tit2) {
                    const Triangle& triangle2 = *tit2;
                #else
                for (int i2=tit1-triangles.begin(); i2<static_cast<int>(triangles.size()); ++i2) {
                    const Triangle& triangle2 = *(triangles.begin()+i2);
                #endif
                    e.Run([&](){
                        matrix(triangle1.index(),triangle2.index()) = base::integrator.integrate(Sfunc,triangle2)*coeff;
                    });
                }
                e.Rethrow();
            }
        }

        template <typename T>
        void N(const double coeff,T& matrix) const {
            if (S_block_is_computed()) {
                N(coeff/Scoeff,matrix,matrix);
            } else {
                SymBloc Sbloc(mesh.triangles().front().index(),mesh.triangles().size());
                S(1.0,Sbloc);
                N(coeff,Sbloc,matrix);
            }
        }

        template <typename T>
        void D(const double coeff,T& matrix) const {
            base::message("D",mesh,mesh);
            base::D(mesh.triangles(),mesh.triangles(),coeff,matrix);
        }

        template <typename T>
        void Dstar(const double coeff,T& matrix) const {
            base::message("D*",mesh,mesh);
            base::D(mesh.triangles(),mesh.triangles(),coeff,matrix);
        }


    private:

        bool S_block_is_computed() const { return Scoeff!=0.0; }

        static double Id(const Triangle& T,const Vertex& V) {
            return (T.contains(V)) ? T.area()/3.0 : 0.0;
        }

        template <typename T1,typename T2>
        void N(const double coeff,const T1& S,T2& matrix) const {

            base::message("N",mesh,mesh);

            ThreadException e;
            ProgressBar pb(mesh.vertices().size());

            // When meshes are equal, optimized computation for a symmetric matrix.

            for (auto vit1=mesh.vertices().begin(); vit1!=mesh.vertices().end(); ++vit1) {
                #pragma omp parallel for
                #if defined NO_OPENMP || defined OPENMP_ITERATOR
                for (auto vit2=vit1; vit2<mesh.vertices().end(); ++vit2) {
                #else
                for (int i2=0;i2<=vit1-mesh.vertices().begin();++i2) {
                    const auto vit2 = mesh.vertices().begin()+i2;
                #endif
                    e.Run([&](){
                        matrix((*vit1)->index(),(*vit2)->index()) += base::N(**vit1,**vit2,mesh,S)*coeff;
                    });
                }
                e.Rethrow();
                ++pb;
            }
        }

        const Mesh&  mesh;
              double Scoeff = 0.0;
    };

    class PartialBlock {
    public:

        PartialBlock(const Mesh& m): mesh(m) { }

        void addD(const double coeff,const Vertices& points,Matrix& matrix) const {
            log_stream(INFORMATION) << "PARTAL OPERATOR D..." << std::endl;
            for (const auto& triangle : mesh.triangles()) {
                const analyticD3 analyD(triangle);
                for (const auto& vertex : points) {
                    const Vect3& integrals = analyD.f(vertex);
                    for (unsigned i=0;i<3;++i)
                        matrix(vertex.index(),triangle.vertex(i).index()) += integrals(i)*coeff;
                }
            }
        }

        void S(const double coeff,const Vertices& points,Matrix& matrix) const {
            log_stream(INFORMATION) << "PARTIAL OPERATOR S..." << std::endl;
            for (const auto& triangle : mesh.triangles()) {
                const analyticS analyS(triangle);
                for (const auto& vertex : points)
                    matrix(vertex.index(),triangle.index()) = coeff*analyS.f(vertex);
            }
        }

    private:

        const Mesh& mesh;
    };

    class NonDiagonalBlock: public BlocksBase  {

        typedef BlocksBase base;

        class Bloc: public Matrix {

            typedef Matrix base;

        public:

            Bloc(const unsigned r0,const unsigned c0,const unsigned n,const unsigned m): base(n,m),i0(r0),j0(c0) { }

            double& operator()(const unsigned i,const unsigned j)       { return base::operator()(i-i0,j-j0); }
            double  operator()(const unsigned i,const unsigned j) const { return base::operator()(i-i0,j-j0); }

        private:

            const unsigned i0;
            const unsigned j0;
        };

    public:

        // This constructor takes the following arguments:
        //  - The 2 interacting meshes.
        //  - The gauss order parameter (for adaptive integration).
        //  - A verbosity parameters (for printing the action on the terminal).

        NonDiagonalBlock(const Mesh& m1,const Mesh& m2,const Integrator& intg): base(intg),mesh1(m1),mesh2(m2) { }

        template <typename T>
        void set_S_block(const double coeff,T& matrix) {
            if (!mesh1.current_barrier() && !mesh2.current_barrier()) {
                S(coeff,matrix);
                Scoeff = coeff;
            }
        }

        template <typename T>
        void set_N_block(const double coeff,T& matrix) const { N(coeff,matrix); }

        template <typename T>
        void set_D_block(const double coeff,T& matrix) const {
            if (!mesh1.current_barrier())
                D(coeff,matrix);
        }

        template <typename T>
        void set_Dstar_block(const double coeff,T& matrix) const {
            if (mesh1!=mesh2 && !mesh2.current_barrier())
                Dstar(coeff,matrix);
        }

        // Various operators take the following arguments:
        //  - The coefficient to be applied to each matrix element (depending on conductivities, ...)
        //  - The storage Matrix for the result

        template <typename T>
        void S(const double coeff,T& matrix) const {
            base::message("S",mesh1,mesh2);
            ThreadException e;
            ProgressBar pb(mesh1.triangles().size());

            // Operator S is given by Sij=\Int G*PSI(I,i)*Psi(J,j) with PSI(l,t) a P0 test function on layer l and triangle t.

            // TODO check the symmetry of S.
            // if we invert tit1 with tit2: results in HeadMat differs at 4.e-5 which is too big.
            // using ADAPT_LHS with tolerance at 0.000005 (for S) drops this at 6.e-6 (but increase the computation time).

            for (const auto& triangle1 : mesh1.triangles()) {

                const analyticS analyS(triangle1);
                const auto& Sfunc = [&analyS](const Vect3& r) { return analyS.f(r); };

                const Triangles& m2_triangles = mesh2.triangles();
                #pragma omp parallel for
                #if defined NO_OPENMP || defined OPENMP_RANGEFOR
                for (const auto& triangle2 : m2_triangles) {
                #elif defined OPENMP_ITERATOR
                for (Triangles::const_iterator tit2=m2_triangles.begin();tit2<m2_triangles.end();++tit2) {
                    const Triangle& triangle2 = *tit2;
                #else
                for (int i2=0;i2<static_cast<int>(m2_triangles.size());++i2) {
                    const Triangle& triangle2 = *(m2_triangles.begin()+i2);
                #endif
                    e.Run([&](){
                        matrix(triangle1.index(),triangle2.index()) = base::integrator.integrate(Sfunc,triangle2)*coeff;
                    });
                }
                e.Rethrow();
                ++pb;
            }
        }

        template <typename T>
        void N(const double coeff,T& matrix) const {
            if (S_block_is_computed()) {
                N(coeff/Scoeff,matrix,matrix);
            } else {
                Bloc Sbloc(mesh1.triangles().front().index(),mesh2.triangles().front().index(),mesh1.triangles().size(),mesh2.triangles().size());
                S(1.0,Sbloc);
                N(coeff,Sbloc,matrix);
            }
        }

        template <typename T>
        void D(const double coeff,T& matrix) const {
            base::message("D",mesh1,mesh2);
            base::D(mesh1.triangles(),mesh2.triangles(),coeff,matrix);
        }

        template <typename T>
        void Dstar(const double coeff,T& matrix) const {
            base::message("D*",mesh1,mesh2);
            base::D(mesh2.triangles(),mesh1.triangles(),coeff,matrix);
        }

    private:

        bool S_block_is_computed() const { return Scoeff!=0.0; }

        template <typename T1,typename T2>
        void N(const double coeff,const T1& S,T2& matrix) const {

            base::message("N",mesh1,mesh2);

            ThreadException e;
            ProgressBar pb(mesh1.vertices().size());
            const VerticesRefs& m2_vertices = mesh2.vertices();
            for (const auto& vertex1 : mesh1.vertices()) {
                #pragma omp parallel for
                #if defined NO_OPENMP || defined OPENMP_RANGEFOR
                for (const auto& vertex2 : m2_vertices) {
                #elif defined OPENMP_ITERATOR
                for (auto vit2=m2_vertices.begin(); vit2<m2_vertices.end(); ++vit2) {
                    const Vertex* vertex2 = *vit2;
                #else
                for (int i2=0; i2<static_cast<int>(m2_vertices.size()); ++i2) {
                    const Vertex* vertex2 = *(m2_vertices.begin()+i2);
                #endif
                    e.Run([&](){
                        matrix(vertex1->index(),vertex2->index()) += base::N(*vertex1,*vertex2,mesh1,mesh2,S)*coeff;
                    });
                }
                e.Rethrow();
                ++pb;
            }
        }

        const Mesh&  mesh1;
        const Mesh&  mesh2;
              double Scoeff = 0.0;
    };

    template <typename BlockType>
    class HeadMatrixBlocks {
    public:

        HeadMatrixBlocks(const BlockType& blk): block(blk) { }

        // SymMatrix is initialized at once, and there is nothing to for blockwise matrix.

        static void init(SymMatrix& matrix) { matrix.set(0.0); }
        void set_blocks(SymMatrix&) const   { }

        // SymmetricBlockMatrix is initialized blockwise.

        static void init(maths::SymmetricBlockMatrix&) { }
        #if 0
        void set_blocks(maths::SymmetricBlockMatrix& matrix) const { // Integrate this in blocks.
            const Ranges& vrange1 = mesh1.vertices_ranges();
            const Ranges& vrange2 = mesh2.vertices_ranges();
            const Ranges& trange1 = { mesh1.triangles_range() };
            const Ranges& trange2 = { mesh2.triangles_range() };

            matrix.add_blocks(trange1,trange2); // S blocks.
            matrix.add_blocks(vrange1,vrange2); // N blocks.
            matrix.add_blocks(trange1,vrange2); // D blocks.
            if (!same_mesh)
                matrix.add_blocks(trange2,vrange1); // D* blocks when they are not the transpose of D blocks.
        }
        #endif

        template <typename T>
        void set_blocks(const double coeffs[3],T& matrix) {
            const double SCondCoeff = coeffs[0];
            const double NCondCoeff = coeffs[1];
            const double DCondCoeff = coeffs[2];
            block.set_S_block(SCondCoeff,matrix);
            block.set_N_block(NCondCoeff,matrix);
            block.set_D_block(DCondCoeff,matrix);
            block.set_Dstar_block(DCondCoeff,matrix);
        }

    private:

        BlockType block;
    };
}
