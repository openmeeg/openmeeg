// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#ifdef WIN32
#define _USE_MATH_DEFINES
#endif

#include <om_common.h>
#include <logger.h>
#include <matrix.h>
#include <symmatrix.h>
#include <matop.h>
#include <geometry.h>
#include <operators.h>
#include <assemble.h>
#include <HeadMatrix.h>

#include <constants.h>

namespace OpenMEEG {

    namespace Details {

        struct AllBlocks {
            AllBlocks() { }
            bool operator()(const Mesh&) const { return false; }
        };

        struct AllButBlock {
            AllButBlock(const Mesh& m): mesh(m) { }
            bool operator()(const Mesh& mesh1) const { return mesh1==mesh; }
            const Mesh& mesh;
        };

        template <typename Selector=AllBlocks()>
        HeadMatrix make_head_matrix(const Geometry& geo,const Integrator& integrator,const Selector& disableBlock=Selector()) {

            log_stream(INFORMATION) << "Assembling Head Matrix" << std::endl;

            HeadMatrix headmat(geo);
            for (const Mesh* mp : geo.interacting_meshes()) {
                const Mesh& mesh = *mp;

                log_stream(INFORMATION) << "Assembling " << mesh.name() << " x " << mesh.name();
                if (disableBlock(mesh)) {
                    log_stream(INFORMATION) << " (skipped)" << std::endl;
                    continue;
                }
                log_stream(INFORMATION) << " (self) setting blocks" << std::endl;

                // TODO: For the scalp the S factor is infinite (one conductivity being 0).
                // What happens in this case ?

                headmat.add_block(mesh,geo.domains(mesh),integrator);
            }

            // Iterate over pairs of adjacent meshes (sharing a domain) to fill the
            // lower half of the HeadMat (since it is symmetric).

            for (const auto& mp : geo.adjacent_meshes()) {
                const Mesh& mesh1 = mp(0);
                const Mesh& mesh2 = mp(1);

                log_stream(INFORMATION) << "Assembling " << mesh1.name() << " x " << mesh2.name() << " setting blocks" << std::endl;

                headmat.add_block(mp,integrator);
            }

            #if 0
            // Deflate all current barriers as one
            log_stream(INFORMATION) << "Deflating current barriers" << std::endl;
            deflate(symmatrix,geo);
            log_stream(INFORMATION) << "done" << std::endl;
            #endif
            return headmat;
        }
    }

    HeadMatrix BlockHeadMat(const Geometry& geo,const Integrator& integrator) {
        return Details::make_head_matrix(geo,integrator,Details::AllBlocks());
    }

    SymMatrix HeadMat(const Geometry& geo,const Integrator& integrator) {

        const HeadMatrix& block_matrix = BlockHeadMat(geo,integrator);
        SymMatrix         head_mat     = block_matrix.make_head_matrix();

        // Deflate all current barriers as one

        log_stream(INFORMATION) << "Deflating current barriers" << std::endl;
        HeadMatrix::deflate(geo,head_mat);
        log_stream(INFORMATION) << "done" << std::endl;

        return head_mat;
    }

    // Construct the head matrix while ignoring the unknown related to the domain cortex (usually the inner domain).

    Matrix PartialHeadMat(const Geometry& geo,const Mesh& cortex,const Integrator& integrator) {

        log_stream(INFORMATION) << "Computing HeadMatrix." << std::endl;
        log_stream(INFORMATION) << "    Found cortex mesh " << cortex.name() << std::endl;
        const HeadMatrix& block_matrix = Details::make_head_matrix(geo,integrator,Details::AllButBlock(cortex));
        const SymMatrix   head_matrix  = block_matrix.make_head_matrix();

        // Copy head_matrix into the returned matrix except for the lines related to the cortex
        // (vertices [i_vb_c, i_ve_c] and triangles [i_tb_c, i_te_c]).

        log_stream(INFORMATION) << "    Copying head matrix into the returned matrix";
        const unsigned Nl = geo.num_unknowns()-cortex.num_unknowns();
        log_stream(INFORMATION) << " with " << Nl << " lines" << std::endl;

        Matrix matrix = Matrix(Nl,head_matrix.ncol());
        matrix = 0.0;
        unsigned iNl = 0;
        for (const auto& mesh : geo.meshes())
            if (mesh!=cortex) {
                log_stream(INFORMATION) << "Processing mesh " << mesh.name() << " : vertices";
                for (const auto& vertex : mesh.vertices())
                    matrix.row(iNl++) = head_matrix.row(vertex->index());
                if (!mesh.is_current_barrier()) {
                    log_stream(INFORMATION) << " : triangles";
                    for (const auto& triangle : mesh.triangles())
                        matrix.row(iNl++) = head_matrix.row(triangle.index());
                }
                log_stream(INFORMATION) << " done." << std::endl;
            }

        return matrix;
    }

    Matrix
    CorticalMat(const Geometry& geo,const SparseMatrix& M,const Matrix& P,const double alpha,const double beta) {
        // Following the article: M. Clerc, J. Kybic "Cortical mapping by Laplace–Cauchy transmission using
        // a boundary element method".
        // M is the head to EEG matrix, while P is a projector

        // Get the gradient of P1&P0 elements on the meshes.

        const unsigned Nc = geo.num_unknowns();
        SymMatrix RR(Nc,Nc);
        RR = 0.0;
        for (const auto& mesh : geo.meshes())
            mesh.gradient_norm2(RR);

        /// Choose Regularization parameter

        const Matrix MM(M.transpose()*M);
        SparseMatrix alphas(Nc,Nc); // diagonal matrix
        Matrix Z;
        double alpha1 = alpha;
        double beta1  = beta;
        if (alpha1<0) { // try an automatic method... TODO find better estimation
            const double nRR_v = RR.submat(0,geo.vertices().size(),0,geo.vertices().size()).frobenius_norm();
            alpha1 = MM.frobenius_norm()/(1.e3*nRR_v);
            beta1  = alpha1*50000.;
            log_stream(INFORMATION) << "AUTOMATIC alphas = " << alpha1 << "\tbeta = " << beta1 << std::endl;
        } else {
            log_stream(INFORMATION) << "alphas = " << alpha << "\tbeta = " << beta << std::endl;
        }

        for (const auto& vertex : geo.vertices())
            alphas(vertex.index(),vertex.index()) = alpha1;

        for (const auto& mesh : geo.meshes())
            if (!mesh.is_current_barrier() )
                for (const auto& triangle : mesh.triangles())
                    alphas(triangle.index(),triangle.index()) = beta1;

        Z = P.transpose()*(MM+alphas*RR)*P;

        // ** PseudoInverse and return **
        // X = P * { (M*P)' * (M*P) + (R*P)' * (R*P) }¡(-1) * (M*P)'m
        // X = P * { P'*M'*M*P + P'*R'*R*P }¡(-1) * P'*M'm
        // X = P * { P'*(MM + a*RR)*P }¡(-1) * P'*M'm
        // X = P * Z¡(-1) * P' * M'm

        const Matrix& rhs = P.transpose()*M.transpose();
        return P*Z.pinverse()*rhs;
    }

    Matrix
    CorticalMat(const Geometry& geo,const SparseMatrix& M,const Matrix& H,const double gamma) {

        // Re-writting of the optimization problem in M. Clerc, J. Kybic "Cortical mapping by Laplace–Cauchy
        // transmission using a boundary element method".
        // with a Lagrangian formulation as in see http://www.math.uh.edu/~rohop/fall_06/Chapter3.pdf eq3.3
        // find argmin(norm(gradient(X)) under constraints:
        // H*X = 0 and M*X = m
        // let G be the gradient norm matrix, l1, l2 the lagrange parameters
        //
        // [ G  H' M'] [ X  ]   [ 0 ]
        // | H  0    | | l1 | = | 0 |
        // [ M     0 ] [ l2 ]   [ m ]
        //
        // {----,----}
        //      K
        // we want a submatrix of the inverse of K (using block-wise inversion, (TODO maybe iterative solution better ?)).
        // Assumptions:
        // - domain_name: the domain containing the sources is an innermost domain (defined as the interior of only one
        //   interface (called Cortex)
        // - Cortex interface is composed of one mesh only (no shared vertices)

        // Concat M to H in HM.

        const unsigned Nl = H.nlin();
        const unsigned Nc = H.ncol();
        const unsigned Nt = Nl+M.nlin();

        Matrix HM(Nt,Nc);
        
        for (unsigned i=0; i<Nl; ++i)
            for (unsigned j=0; j<Nc; ++j)
                HM(i,j) = H(i,j);

        for (unsigned i=Nl;i<Nt;++i)
            for (unsigned j=0;j<Nc;++j)
                HM(i,j) = M(i-Nl,j);

        // Get the gradient of P1&P0 elements on the meshes.

        SymMatrix G(Nc);
        G = 0.0;
        for (const auto& mesh : geo.meshes())
            mesh.gradient_norm2(G);

        // Multiply by gamma the submat of current gradient norm2

        for (const auto& mesh : geo.meshes())
            if (!mesh.is_current_barrier())
                for (const auto& triangle1 : mesh.triangles())
                    for (const auto& triangle2 : mesh.triangles())
                        G(triangle1.index(),triangle2.index()) *= gamma;

        log_stream(INFORMATION) << "gamma = " << gamma << std::endl;

        G.invert();
        return (G*HM.transpose()*(HM*G*HM.transpose()).inverse()).submat(0,Nc,Nl,M.nlin());
    }

    Matrix Surf2VolMat(const Geometry& geo,const Matrix& points) {

        // Find the points per domain and generate the indices for the m_points
        // What happens if a point is on the boundary of a domain ? TODO

        std::map<const Domain*,Vertices> m_points;
        unsigned index = 0;
        for (unsigned i=0; i<points.nlin(); ++i) {
            const Vect3 point(points(i,0),points(i,1),points(i,2));
            const Domain& domain = geo.domain(point);
            if (domain.conductivity()==0.0) {
                log_stream(INFORMATION) << " Surf2Vol: Point [ " << points.row(i) << "]"
                                        << " is inside a non-conductive domain. Point is dropped." << std::endl;
            } else {
                m_points[&domain].push_back(Vertex(point,index++));
            }
        }

        // At this point index is the total number of inside points.

        const unsigned size = index; // total number of inside points

        Matrix mat(size,geo.num_unknowns());
        mat = 0.0;

        for (const auto& [domainptr,pts] : m_points) {
            for (const auto& boundary : domainptr->boundaries())
                for (const auto& omesh : boundary.interface().oriented_meshes()) {
                    const Mesh& mesh = omesh.mesh();
                    const PartialBlockWriter writer(mesh);
                    const double coeff = boundary.mesh_orientation(omesh)*K;
                    writer.addD(-coeff,pts,mat);
                    if (!mesh.is_current_barrier())
                        writer.S(coeff/domainptr->conductivity(),pts,mat);
                }
        }

        return mat;
    }
}
