// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <map>
#include <functional>

#include <matrix.h>
#include <symmatrix.h>
#include <block.h>
#include <range.h>
#include <mesh.h>
#include <ConductivityCoefficients.h>
#include <operators.h>

namespace OpenMEEG {

    namespace Details {

        // This class describes how a specific mesh matrix is mapped into the
        // full head matrix. Triangles are contiguous in each mesh so their indices
        // only depends on one shift index. Vertices can be shared across different
        // meshes so require a full index map.

        class MeshMapping {
        public:

            MeshMapping(const Mesh& m):
                t_range(m.triangles_range()),
                v_indmap(m.vertices_mapping()),v_indmap_inv(m.inverse_vertices_mapping(v_indmap))
            { }

            size_t num_triangles() const { return t_range.length(); }
            size_t num_vertices()  const { return v_indmap.size();  }

            size_t triangle_mesh_to_head_matrix(const size_t i) const { return i+t_range.start(); }
            size_t triangle_head_matrix_to_mesh(const size_t i) const { return i-t_range.start(); }

            size_t vertex_mesh_to_head_matrix(const size_t i) const { return v_indmap[i];        }
            size_t vertex_head_matrix_to_mesh(const size_t i) const { return v_indmap_inv.at(i); }

        private:

            using IndexMap    = Mesh::VertexMapping;
            using IndexInvMap = std::map<size_t,size_t>;
            using TRange      = OpenMEEG::maths::Range;

            const TRange      t_range;      // Index range for the triangles in the general symmetric matrix.
            const IndexMap    v_indmap;     // Mapping of this mesh vertices in the general symmetric matrix.
            const IndexInvMap v_indmap_inv; // Inverse mapping of v_indmap.
        };

        // Helper classes to fill the global symmetric head matrix from sub-blocks.

        template <typename Matrix> class BlockWriter;

        template <>
        class BlockWriter<const Matrix>: public Block<const Matrix> {

            using base = Block<const Matrix>;

        public:

            using base::base;

            void add_block_contribution(const double coeff,SymMatrix& head_matrix) const {
                for (unsigned i=0; i<matrix.nlin(); ++i)
                    for (unsigned j=0; j<matrix.ncol(); ++j)
                        head_matrix(ind1(i),ind2(j)) += coeff*matrix(i,j);
            }
        };

        template <>
        class BlockWriter<const SymMatrix>: public Block<const SymMatrix> {

            using base = Block<const SymMatrix>;

        public:

            using base::base;

            void add_block_contribution(const double coeff,SymMatrix& head_matrix) const {
                for (unsigned i=0; i<matrix.nlin(); ++i)
                    for (unsigned j=i; j<matrix.ncol(); ++j)
                        head_matrix(ind1(i),ind2(j)) += coeff*matrix(i,j);
            }
        };

        // The next two classes describe the various mesh blocks for the symmetric approach.
        // Each block is of the form
        //   N | D*    Where N is of size V1 x V2, S of size T1 x T2, D of size T1 x V2 and D* of size V1 x T2,
        //   ------    with V1 and V2 being the number of vertices of mesh1 and mesh2 respectively,
        //   D | S     and T1 and T2 being the number of triangles  of mesh1 and mesh2 respectively.

        // When the two meshes are the same (diagonal blocks), N and S are symmetric matrices and D* = D^T,
        // so there is no need to tore it independently.

        class DiagonalBlock {
        public:

            // TODO: Optimize memory by not allocating memory for S and D in cases of current barriers.

            DiagonalBlock(const MeshMapping& map,const DomainsRefs& doms,const bool mesh_is_current_barrier):
                mapping(map),domains(doms),N_block(map.num_vertices())
            {
                N_block = 0.0;
                if (!mesh_is_current_barrier) {
                    S_block.resize(map.num_triangles());
                    D_block.resize(map.num_triangles(),map.num_vertices());
                    D_block = 0.0;
                }
            }

            Block<SymMatrix> N() {
                const auto vmap = [&](const size_t i) { return mapping.vertex_head_matrix_to_mesh(i); };
                return Block<SymMatrix>(N_block,vmap,vmap);
            }

            Block<SymMatrix> S() {
                const auto tmap = [&](const size_t i) { return mapping.triangle_head_matrix_to_mesh(i); };
                return Block<SymMatrix>(S_block,tmap,tmap);
            }

            Block<Matrix> D() {
                const auto vmap = [&](const size_t i) { return mapping.vertex_head_matrix_to_mesh(i);   };
                const auto tmap = [&](const size_t i) { return mapping.triangle_head_matrix_to_mesh(i); };
                return Block<Matrix>(D_block,tmap,vmap);
            }

            void update_block(SymMatrix& head_matrix) const {
                const auto vmap = [&](const size_t i) { return mapping.vertex_mesh_to_head_matrix(i);   };
                const auto tmap = [&](const size_t i) { return mapping.triangle_mesh_to_head_matrix(i); };

                const ConductivityCoefficients cond_coeffs(domains);

                BlockWriter<const SymMatrix> N(N_block,vmap,vmap);
                N.add_block_contribution(cond_coeffs.N(),head_matrix);

                if (!S_block.empty()) {
                    BlockWriter<const SymMatrix> S(S_block,tmap,tmap);
                    S.add_block_contribution(cond_coeffs.S(),head_matrix);
                }

                if (!D_block.empty()) {
                    BlockWriter<const Matrix> D(D_block,tmap,vmap);
                    D.add_block_contribution(cond_coeffs.D(),head_matrix);
                }
            }

        private:

            const MeshMapping& mapping; // Mapping of vertex and triangle indices.
            const DomainsRefs  domains; // Domains that the mesh of this block delimits.

            // Blocks.

            SymMatrix N_block;
            SymMatrix S_block;
            Matrix    D_block;
        };

        class NonDiagonalBlock {
        public:

            NonDiagonalBlock(const MeshMapping& map1,const MeshMapping& map2,const MeshPair& mp):
                mapping1(map1),mapping2(map2),mesh_pair(mp),N_block(map1.num_vertices(),map2.num_vertices())
            {
                const Mesh& mesh1 = mp(0);
                const Mesh& mesh2 = mp(1);

                N_block = 0.0;
                if (!mesh1.is_current_barrier()) {
                    D_block.resize(map1.num_triangles(),map2.num_vertices());
                    D_block = 0.0;
                    if (!mesh2.is_current_barrier())
                        S_block.resize(map1.num_triangles(),map2.num_triangles());
                }
                if (!mesh2.is_current_barrier()) {
                    D_star_block.resize(map1.num_vertices(),map2.num_triangles());
                    D_star_block = 0.0;
                }
            }

            Block<Matrix> N() {
                const auto vmap1 = [&](const size_t i) { return mapping1.vertex_head_matrix_to_mesh(i); };
                const auto vmap2 = [&](const size_t i) { return mapping2.vertex_head_matrix_to_mesh(i); };
                return Block<Matrix>(N_block,vmap1,vmap2);
            }

            Block<Matrix> S() {
                const auto tmap1 = [&](const size_t i) { return mapping1.triangle_head_matrix_to_mesh(i); };
                const auto tmap2 = [&](const size_t i) { return mapping2.triangle_head_matrix_to_mesh(i); };
                return Block<Matrix>(S_block,tmap1,tmap2);
            }

            Block<Matrix> D() {
                const auto tmap = [&](const size_t i) { return mapping1.triangle_head_matrix_to_mesh(i); };
                const auto vmap = [&](const size_t i) { return mapping2.vertex_head_matrix_to_mesh(i); };
                return Block<Matrix>(D_block,tmap,vmap);
            }

            Block<Matrix> Dstar() {
                const auto vmap = [&](const size_t i) { return mapping1.vertex_head_matrix_to_mesh(i); };
                const auto tmap = [&](const size_t i) { return mapping2.triangle_head_matrix_to_mesh(i); };
                return Block<Matrix>(D_star_block,vmap,tmap);
            }

            void update_block(SymMatrix& head_matrix) const {
                const auto vmap1 = [&](const size_t i) { return mapping1.vertex_mesh_to_head_matrix(i);   };
                const auto vmap2 = [&](const size_t i) { return mapping2.vertex_mesh_to_head_matrix(i);   };
                const auto tmap1 = [&](const size_t i) { return mapping1.triangle_mesh_to_head_matrix(i); };
                const auto tmap2 = [&](const size_t i) { return mapping2.triangle_mesh_to_head_matrix(i); };

                const ConductivityCoefficients cond_coeffs(mesh_pair);

                BlockWriter<const Matrix> N(N_block,vmap1,vmap2);
                N.add_block_contribution(cond_coeffs.N(),head_matrix);

                if (!S_block.empty()) {
                    BlockWriter<const Matrix> S(S_block,tmap1,tmap2);
                    S.add_block_contribution(cond_coeffs.S(),head_matrix);
                }

                if (!D_block.empty()) {
                    BlockWriter<const Matrix> D(D_block,tmap1,vmap2);
                    D.add_block_contribution(cond_coeffs.D(),head_matrix);
                }

                if (!D_star_block.empty()) {
                    BlockWriter<const Matrix> D_star(D_star_block,vmap1,tmap2);
                    D_star.add_block_contribution(cond_coeffs.D(),head_matrix);
                }
            }

        private:

            const MeshMapping& mapping1;  // Mapping of vertex and triangle indices for mesh 1.
            const MeshMapping& mapping2;  // Mapping of vertex and triangle indices for mesh 2.
            const MeshPair&    mesh_pair; // Domains that the mesh of this block delimits.

            // Blocks.

            Matrix N_block;
            Matrix S_block;
            Matrix D_block;
            Matrix D_star_block;
        };
    }

    // A representation of the overall head matrix in which each block with unknown conductivity
    // is kept separate. This allows a fast computation of several matrices with different
    // conductivities.

    class HeadMatrix {
    public:

        HeadMatrix(const Geometry& geom): headmat(geom.num_unknowns()) { headmat = 0.0; }

        // Diagonal block.

        void add_block(const Mesh& mesh,const DomainsRefs& doms,const Integrator& integrator) {
            if (has_unknown_conductivity(doms)) {
                using namespace Details;
                const MeshMapping& map1 = add_mapping(mesh);
                DiagonalBlock& block = diag_mesh_blocks.emplace(&mesh,DiagonalBlock(map1,doms,mesh.is_current_barrier())).first->second;
                DiagonalBlockWriter block_writer(mesh,integrator);
                set_blocks(block_writer,block);
            } else {
                DiagonalBlockWriter block(mesh,integrator);
                set_blocks(ConductivityCoefficients(doms),block);
            }
        }

        // Non diagonal block.

        void add_block(const MeshPair& mp,const Integrator& integrator) {
            using namespace Details;

            const Mesh& mesh1 = mp(0);
            const Mesh& mesh2 = mp(1);

            if (has_unknown_conductivity(mp.common_domains())) {
                const MeshMapping& map1 = add_mapping(mesh1);
                const MeshMapping& map2 = add_mapping(mesh2);
                NonDiagonalBlock& block = non_diag_mesh_blocks.emplace(mp,NonDiagonalBlock(map1,map2,mp)).first->second;
                NonDiagonalBlockWriter block_writer(mesh1,mesh2,integrator);
                set_blocks(block_writer,block);
                set_Dstar_block(block_writer,block);
            } else {
                const ConductivityCoefficients cond_coeffs(mp);
                NonDiagonalBlockWriter block(mesh1,mesh2,integrator);
                set_blocks(cond_coeffs,block);
                set_Dstar_block(cond_coeffs.D(),block);
            }
        }

        SymMatrix make_head_matrix() const {
            if (diag_mesh_blocks.empty() && non_diag_mesh_blocks.empty())
                return headmat;

            SymMatrix result = headmat;
            for (const auto& [index,diag_block] : diag_mesh_blocks)
                diag_block.update_block(result);
            for (const auto& [index,non_diag_block] : non_diag_mesh_blocks)
                non_diag_block.update_block(result);
            return result;
        }

        // Deflate all current barriers as one.

        static void deflate(const Geometry& geo,SymMatrix& headmat) {
            ThreadException e;
            for (const auto& part : geo.isolated_parts()) {
                unsigned num_vertices = 0;
                unsigned i_first = 0;
                for (const auto& meshptr : part)
                    if (meshptr->is_external()){
                        num_vertices += meshptr->vertices().size();
                        if (i_first==0)
                            i_first = meshptr->vertices().front()->index();
                    }
                const double coef = headmat(i_first,i_first)/num_vertices;
                for (const auto& meshptr : part)
                    if (meshptr->is_external()) {
                        const auto& vertices = meshptr->vertices();
                        for (auto vit1=vertices.begin(); vit1!=vertices.end(); ++vit1) {
                            #pragma omp parallel for
                            #if defined NO_OPENMP || defined OPENMP_ITERATOR
                            for (auto vit2=vit1; vit2<vertices.end(); ++vit2) {
                            #else
                            for (int i2=vit1-vertices.begin();i2<static_cast<int>(vertices.size());++i2) {
                                const auto vit2 = vertices.begin()+i2;
                            #endif
                                e.run([&]() {
                                    headmat((*vit1)->index(),(*vit2)->index()) += coef;
                                });
                            }
                            e.rethrow();
                        }
                    }
            }
        }

    private:

        static bool has_unknown_conductivity(const DomainsRefs& doms) {
            for (const Domain* domain : doms)
                if (!domain->has_conductivity())
                    return true;
            return false;
        }

        Details::MeshMapping& add_mapping(const Mesh& mesh) {
            MeshMappings::iterator it = mesh_mappings.find(&mesh);
            if (it!=mesh_mappings.end())
                return it->second;
            return mesh_mappings.emplace(&mesh,Details::MeshMapping(mesh)).first->second;
        }

        template <typename BlockWriter,typename Block>
        void set_blocks(BlockWriter& writer,Block& matrix) {
            auto S = matrix.S();
            if (!S.empty())
                writer.set_S_block(1.0,S);

            auto N = matrix.N();
            writer.set_N_block(1.0,S,N);

            auto D = matrix.D();
            if (!D.empty())
                writer.set_D_block(1.0,D);
        }

        void set_Dstar_block(NonDiagonalBlockWriter& block,Details::NonDiagonalBlock& matrix) {
            auto Dstar = matrix.Dstar().transpose();
            if (!Dstar.empty())
                block.set_Dstar_block(1.0,Dstar);
        }

        template <typename BlockWriter>
        void set_blocks(const ConductivityCoefficients& cond_coeffs,BlockWriter& writer) {
            writer.set_S_block(cond_coeffs.S(),headmat);
            writer.set_N_block(cond_coeffs.N(),headmat);
            writer.set_D_block(cond_coeffs.D(),headmat);
        }

        void set_Dstar_block(const double coeff,NonDiagonalBlockWriter& writer) {
            writer.set_Dstar_block(coeff,headmat);
        }

        using MeshIndex         = const Mesh*;
        using MeshMappings      = std::map<MeshIndex,Details::MeshMapping>;
        using DiagonalBlocks    = std::map<MeshIndex,Details::DiagonalBlock>;
        using NonDiagonalBlocks = std::map<MeshPair,Details::NonDiagonalBlock>;

        MeshMappings      mesh_mappings;        // Vertex/triangle mappings that allow to place individual block in the full head matrix.
        DiagonalBlocks    diag_mesh_blocks;     // Diagonal blocks involving a mesh bordering a domain with unknown conductivity.
        NonDiagonalBlocks non_diag_mesh_blocks; // Non diagonal blocks involving a mesh bordering a domain with unknown conductivity.
        SymMatrix         headmat;              // The full matrix contains everything that does not depend on a variable conductivity.
    };
}
