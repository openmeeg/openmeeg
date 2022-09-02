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
#include <matrix.h>
#include <symmatrix.h>
#include <matop.h>
#include <geometry.h>
#include <operators.h>
#include <assemble.h>

#include <constants.h>

namespace OpenMEEG {

    namespace Details {

        template <typename T>
        void deflate(T& M,const Geometry& geo) {
            //  deflate all current barriers as one
            for (const auto& part : geo.isolated_parts()) {
                unsigned nb_vertices = 0;
                unsigned i_first = 0;
                for (const auto& meshptr : part)
                    if (meshptr->outermost()){
                        nb_vertices += meshptr->vertices().size();
                        if (i_first==0)
                            i_first = meshptr->vertices().front()->index();
                    }
                const double coef = M(i_first,i_first)/nb_vertices;
                for (const auto& meshptr : part)
                    if (meshptr->outermost()) {
                        const auto& vertices = meshptr->vertices();
                        for (auto vit1=vertices.begin(); vit1!=vertices.end(); ++vit1) {
                            #pragma omp parallel for
                            #if defined NO_OPENMP || defined OPENMP_ITERATOR
                            for (auto vit2=vit1; vit2<vertices.end(); ++vit2) {
                            #else
                            for (int i2=vit1-vertices.begin();i2<static_cast<int>(vertices.size());++i2) {
                                const auto vit2 = vertices.begin()+i2;
                            #endif
                                M((*vit1)->index(),(*vit2)->index()) += coef;
                            }
                        }
                    }
            }
        }

        struct AllBlocks {
            AllBlocks() { }
            bool operator()(const Mesh&,const Mesh&) const { return false; }
        };

        struct AllButBlock {
            AllButBlock(const Mesh& m): mesh(m) { }
            bool operator()(const Mesh& mesh1,const Mesh& mesh2) const { return mesh1==mesh2 && mesh1==mesh; }
            const Mesh& mesh;
        };

        template <typename TYPE,typename Selector>
        TYPE HeadMatrix(const Geometry& geo,const Integrator& integrator,const Selector& disableBlock,const bool verbose=true) {

            TYPE symmatrix(geo.nb_parameters()-geo.nb_current_barrier_triangles());
            HeadMatrixBlocks<TYPE>::init(symmatrix);

            // Iterate over pairs of communicating meshes (sharing a domains) to fill the
            // lower half of the HeadMat (since it is symmetric).

            for (const auto& mp : geo.communicating_mesh_pairs()) {
                const Mesh& mesh1 = mp(0);
                const Mesh& mesh2 = mp(1);

                if (disableBlock(mesh1,mesh2))
                    continue;

                const double factor     =  mp.relative_orientation()*K;
                const double SCondCoeff =  factor*geo.sigma_inv(mesh1,mesh2);
                const double NCondCoeff =  factor*geo.sigma(mesh1,mesh2);
                const double DCondCoeff = -factor*geo.indicator(mesh1,mesh2);

                const double coeffs[3] = { SCondCoeff, NCondCoeff, DCondCoeff };

                if (&mesh1==&mesh2) {
                    HeadMatrixBlocks<DiagonalBlock> operators(DiagonalBlock(mesh1,integrator,verbose));
                    operators.set_blocks(coeffs,symmatrix);
                } else {
                    HeadMatrixBlocks<NonDiagonalBlock> operators(NonDiagonalBlock(mesh1,mesh2,integrator,verbose));
                    operators.set_blocks(coeffs,symmatrix);
                }
            }

            // Deflate all current barriers as one

            deflate(symmatrix,geo);

            return symmatrix;
        }
    }

    #if 0
    SymMatrix conductivity_coefficients(const Geometry& geo) {
        SymMatrix cond_coeffs(10,10);
        for (const auto& mp : geo.communicating_mesh_pairs()) {
            const Mesh& mesh1 = mp(0);
            const Mesh& mesh2 = mp(1);
            const double factor = mp.relative_orientation()*K;
            const double Scoeff =  factor*geo.sigma_inv(mesh1,mesh2);
            const double Ncoeff =  factor*geo.sigma(mesh1,mesh2);
            const double Dcoeff = -factor*geo.indicator(mesh1,mesh2);
        }
        return cond_coeffs;
    }
    #endif

    SymMatrix HeadMat(const Geometry& geo,const Integrator& integrator,const bool verbose) {
        return Details::HeadMatrix<SymMatrix>(geo,integrator,Details::AllBlocks(),verbose);
    }

    Matrix HeadMatrix(const Geometry& geo,const Interface& Cortex,const Integrator& integrator,const unsigned extension=0,const bool verbose=true) {

        const Mesh& cortex = Cortex.oriented_meshes().front().mesh();
        const SymMatrix& symmatrix = Details::HeadMatrix<SymMatrix>(geo,integrator,Details::AllButBlock(cortex),verbose);

        // Copy symmatrix into the returned matrix except for the lines related to the cortex
        // (vertices [i_vb_c, i_ve_c] and triangles [i_tb_c, i_te_c]).

        const unsigned Nl = geo.nb_parameters()-geo.nb_current_barrier_triangles()-Cortex.nb_vertices()-Cortex.nb_triangles()+extension;

        Matrix matrix = Matrix(Nl,symmatrix.ncol());
        matrix.set(0.0);
        unsigned iNl = 0;
        for (const auto& mesh : geo.meshes())
            if (mesh!=cortex) {
                for (const auto& vertex : mesh.vertices())
                    matrix.setlin(iNl++,symmatrix.getlin(vertex->index()));
                if (!mesh.current_barrier())
                    for (const auto& triangle : mesh.triangles())
                        matrix.setlin(iNl++,symmatrix.getlin(triangle.index()));
            }

        return matrix;
    }

    Matrix CorticalMat(const Geometry& geo,const SparseMatrix& M,const std::string& domain_name,
                       const double alpha,const double beta,const std::string& filename,const Integrator& integrator,
                       const bool verbose)
    {
        // Following the article: M. Clerc, J. Kybic "Cortical mapping by Laplace–Cauchy transmission using
        // a boundary element method".
        // Assumptions:
        // - domain_name: the domain containing the sources is an innermost domain (defined as the interior of
        //   only one interface (called Cortex))
        // - Cortex interface is composed of one mesh only (no shared vertices).

        const Domain& SourceDomain = geo.domain(domain_name);
        const Interface& Cortex    = SourceDomain.boundaries().front().interface();

        om_error(SourceDomain.boundaries().size()==1);
        om_error(Cortex.oriented_meshes().size()==1);

        // shape of the new matrix:

        Matrix P;
        std::fstream f(filename.c_str());
        if (!f) {
            const Matrix& mat = HeadMatrix(geo,Cortex,integrator,verbose);

            //  Construct P: the null-space projector.
            //  P is a projector: P^2 = P and mat*P*X = 0

            P = nullspace_projector(mat);
            if (filename.length()!=0) {
                if (verbose)
                    std::cout << "Saving projector P (" << filename << ")." << std::endl;
                P.save(filename);
            }
        } else {
            if (verbose)
                std::cout << "Loading projector P (" << filename << ")." << std::endl;
            P.load(filename);
        }

        // ** Get the gradient of P1&P0 elements on the meshes **

        const unsigned Nc = geo.nb_parameters()-geo.nb_current_barrier_triangles();
        SymMatrix RR(Nc,Nc);
        RR.set(0.0);
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
            alphas.set(0.);
            alpha1 = MM.frobenius_norm()/(1.e3*nRR_v);
            beta1  = alpha1*50000.;
            if (verbose)
                std::cout << "AUTOMATIC alphas = " << alpha1 << "\tbeta = " << beta1 << std::endl;
        } else {
            if (verbose)
                std::cout << "alphas = " << alpha << "\tbeta = " << beta << std::endl;
        }

        for (const auto& vertex : geo.vertices())
            alphas(vertex.index(),vertex.index()) = alpha1;

        for (const auto& mesh : geo.meshes())
            if (!mesh.current_barrier() )
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

    Matrix CorticalMat2(const Geometry& geo,const SparseMatrix& M,const std::string& domain_name,
                        const double gamma,const std::string& filename,const Integrator& integrator,
                        const bool verbose)
    {
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

        const Domain& SourceDomain = geo.domain(domain_name);
        const Interface& Cortex    = SourceDomain.boundaries().front().interface();

        om_error(SourceDomain.boundaries().size()==1);
        om_error(Cortex.oriented_meshes().size()==1);

        // shape of the new matrix:

        std::fstream f(filename.c_str());
        Matrix H;
        if (!f) {
            H = HeadMatrix(geo,Cortex,integrator,M.nlin(),verbose);
            if (filename.length()!=0) {
                if (verbose)
                    std::cout << "Saving matrix H (" << filename << ")." << std::endl;
                H.save(filename);
            }
        } else {
            if (verbose)
                std::cout << "Loading matrix H (" << filename << ")." << std::endl;
            H.load(filename);
        }

        // Why don't we save the full matrix H instead of only the upper block ??
        // Concat M to H
        // TODO TODO TODO: Integrate in HeadMatrix...

        const unsigned Nl = H.nlin()-M.nlin();
        const unsigned Nc = H.ncol();
        for (unsigned i=Nl;i<H.nlin();++i)
            for (unsigned j=0;j<Nc;++j)
                H(i,j) = M(i-Nl,j);

        // Get the gradient of P1&P0 elements on the meshes.

        SymMatrix G(Nc);
        G.set(0.0);
        for (const auto& mesh : geo.meshes())
            mesh.gradient_norm2(G);

        // Multiply by gamma the submat of current gradient norm2

        for (const auto& mesh : geo.meshes())
            if (!mesh.current_barrier())
                for (const auto& triangle1 : mesh.triangles())
                    for (const auto& triangle2 : mesh.triangles())
                        G(triangle1.index(),triangle2.index()) *= gamma;

        if (verbose)
            std::cout << "gamma = " << gamma << std::endl;

        G.invert();
        return (G*H.transpose()*(H*G*H.transpose()).inverse()).submat(0,Nc,Nl,M.nlin());
    }

    Matrix Surf2VolMat(const Geometry& geo,const Matrix& points,const bool verbose) {

        // Find the points per domain and generate the indices for the m_points
        // What happens if a point is on the boundary of a domain ? TODO

        std::map<const Domain*,Vertices> m_points;
        unsigned index = 0;
        for (unsigned i=0; i<points.nlin(); ++i) {
            const Vect3 point(points(i,0),points(i,1),points(i,2));
            const Domain& domain = geo.domain(point);
            if (domain.conductivity()==0.0) {
                std::cerr << " Surf2Vol: Point [ " << points.getlin(i) << "]"
                          << " is inside a non-conductive domain. Point is dropped." << std::endl;
            } else {
                m_points[&domain].push_back(Vertex(point,index++));
            }
        }

        // At this point index is the total number of inside points.

        const unsigned size = index; // total number of inside points

        Matrix mat(size,geo.nb_parameters()-geo.nb_current_barrier_triangles());
        mat.set(0.0);

        for (const auto& [domainptr,pts] : m_points) {
            for (const auto& boundary : domainptr->boundaries())
                for (const auto& omesh : boundary.interface().oriented_meshes()) {
                    const Mesh& mesh = omesh.mesh();
                    const PartialBlock block(mesh, verbose);
                    const double coeff = boundary.mesh_orientation(omesh)*K;
                    block.addD(-coeff,pts,mat);
                    if (!mesh.current_barrier())
                        block.S(coeff/domainptr->conductivity(),pts,mat);
                }
        }

        return mat;
    }
}
