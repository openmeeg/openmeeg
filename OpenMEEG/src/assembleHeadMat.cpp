/* OpenMEEG

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

#if WIN32
#define _USE_MATH_DEFINES
#endif

#include <om_common.h>
#include <matrix.h>
#include <symmatrix.h>
#include <geometry.h>
#include <operators.h>
#include <assemble.h>

namespace OpenMEEG {

    template <typename T>
    void deflate(T& M,const Interface& i,const double coef) {
        //  deflate the Matrix
        for (Interface::const_iterator omit=i.begin();omit!=i.end();++omit) {
            for (auto vit1=omit->mesh().vertices().begin();vit1!=omit->mesh().vertices().end();++vit1) {
                #pragma omp parallel for
                for (auto vit2=vit1;vit2<omit->mesh().vertices().end();++vit2)
                    M((*vit1)->index(),(*vit2)->index()) += coef;
            }
        }
    }

    template <typename T>
    void deflate(T& M,const Geometry& geo) {
        //  deflate all current barriers as one
        unsigned nb_vertices=0,i_first=0;
        double coef=0.0;
        for ( std::vector<Strings>::const_iterator git = geo.geo_group().begin(); git != geo.geo_group().end(); ++git) {
            nb_vertices=0;
            i_first=0;
            for(std::vector<std::string>::const_iterator mit=git->begin();mit!=git->end();++mit){
                const Mesh msh=geo.mesh(*mit);
                if(msh.outermost()){
                    nb_vertices += msh.vertices().size();
                    if(i_first==0)
                        i_first=(*msh.vertices().begin())->index();
                }
            }
            coef=M(i_first,i_first)/nb_vertices;
            for(std::vector<std::string>::const_iterator mit=git->begin();mit!=git->end();++mit){
                Mesh msh=geo.mesh(*mit);
                if(msh.outermost())
                    for(auto vit1=msh.vertices().begin();vit1!=msh.vertices().end();++vit1) {
                        #pragma omp parallel for
                        for (auto vit2=vit1;vit2<msh.vertices().end();++vit2)
                            M((*vit1)->index(),(*vit2)->index()) += coef;
                    }
            }
        }
    }

    HeadMat::HeadMat(const Geometry& geo,const unsigned gauss_order) {

        SymMatrix& mat = *this;

        mat = SymMatrix((geo.nb_parameters()-geo.nb_current_barrier_triangles()));
        mat.set(0.0);
        double K = 1.0 / (4.0 * M_PI);

        // We iterate over the meshes (or pair of domains) to fill the lower half of the HeadMat (since its symmetry)

        for (const auto& mesh1 : geo.meshes())
            if (!mesh1.isolated()) {
                for (const auto& mesh2 : geo.meshes()) {
                    if ((!mesh2.isolated()) && (geo.sigma(mesh1,mesh2)!=0.0)) {
                        // if mesh1 and mesh2 communicate, i.e they are used for the definition of a common domain
                        const int orientation = geo.oriented(mesh1,mesh2); // equals  0, if they don't have any domains in common
                                                                           // equals  1, if they are both oriented toward the same domain
                                                                           // equals -1, if they are not
                        if (orientation!=0) {
                            double Scoeff =  orientation*geo.sigma_inv(mesh1,mesh2)*K;
                            double Dcoeff = -orientation*geo.indicator(mesh1,mesh2)*K;
                            double Ncoeff;

                            if ((!mesh1.current_barrier()) && (!mesh2.current_barrier()) ) {
                                // Computing S block first because it's needed for the corresponding N block
                                operatorS(mesh1,mesh2,mat,Scoeff,gauss_order);
                                Ncoeff = geo.sigma(mesh1,mesh2)/geo.sigma_inv(mesh1,mesh2);
                            }else{
                                Ncoeff = orientation*geo.sigma(mesh1,mesh2)*K;
                            }

                            if (!mesh1.current_barrier()){
                                // Computing D block
                                operatorD(mesh1,mesh2,mat,Dcoeff,gauss_order,false);
                            }
                            if ((mesh1!=mesh2) && (!mesh2.current_barrier())){
                                // Computing D* block
                                operatorD(mesh1,mesh2,mat,Dcoeff,gauss_order,true);
                            }
                            // Computing N block
                            operatorN(mesh1,mesh2,mat,Ncoeff,gauss_order);
                        }
                    }

                    //  Loop only on oriented pairs of meshes.

                    if (&mesh2==&mesh1)
                        break;
                }
            }

        // Deflate all current barriers as one

        deflate(mat,geo);
    }

    CorticalMat::CorticalMat(const Geometry& geo,const Head2EEGMat& M,const std::string& domain_name,
                             const unsigned gauss_order,const double alpha,const double beta,const std::string& filename)
    {
        Matrix& mat = *this;

        // Following the article: M. Clerc, J. Kybic "Cortical mapping by Laplace–Cauchy transmission using a boundary element method".
        // Assumptions:
        // - domain_name: the domain containing the sources is an innermost domain (defined as the interior of only one interface (called Cortex))
        // - Cortex interface is composed of one mesh only (no shared vertices).

        const Domain& SourceDomain = geo.domain(domain_name);
        const Interface& Cortex    = SourceDomain.boundaries().front().interface();
        const Mesh& cortex         = Cortex.front().mesh();
        
        om_error(SourceDomain.boundaries().size()==1);
        om_error(Cortex.size()==1);

        // shape of the new matrix:

        const unsigned Nl = geo.nb_parameters()-geo.nb_current_barrier_triangles()-Cortex.nb_vertices()-Cortex.nb_triangles();
        const unsigned Nc = geo.nb_parameters()-geo.nb_current_barrier_triangles();
        std::fstream f(filename.c_str());
        Matrix P;
        if (!f) {

            // Build the HeadMat:
            // The following is the same as assemble_HM except N_11, D_11 and S_11 are not computed. TODO ?

            SymMatrix mat_temp(Nc);
            mat_temp.set(0.0);
            double K = 1.0 / (4.0 * M_PI);

            // Iterate over the meshes (or pair of domains) to fill the lower half of the HeadMat (since its symmetry)

            for (const auto& mesh1 : geo.meshes())
                for (const auto& mesh2 : geo.meshes()) {

                    // If mesh1 and mesh2 communicate, i.e they are used for the definition of a common domain

                    const int orientation = geo.oriented(mesh1,mesh2);

                    // orientation equals 0, if the 2 meshes don't have any domains in common
                    // Equals  1, if they are both oriented toward the same domain
                    // Equals -1, otherwise.

                    if (orientation!=0) { //    Interacting meshes.
                        const double Scoeff =  orientation*geo.sigma_inv(mesh1,mesh2)*K;
                        const double Dcoeff = -orientation*geo.indicator(mesh1,mesh2)*K;
                        double Ncoeff;
                        if (!(mesh1.current_barrier() || mesh2.current_barrier()) && ((mesh1!=mesh2) || (mesh1!=cortex))) {
                            // Computing S block first because it's needed for the corresponding N block
                            operatorS(mesh1,mesh2,mat_temp,Scoeff,gauss_order);
                            Ncoeff = geo.sigma(mesh1,mesh2)/geo.sigma_inv(mesh1,mesh2);
                        } else {
                            Ncoeff = orientation*geo.sigma(mesh1,mesh2)*K;
                        }

                        if (!mesh1.current_barrier() && (((mesh1!=mesh2) || (mesh1!=cortex)))) // Computing D block
                            operatorD(mesh1,mesh2,mat_temp,Dcoeff,gauss_order,false);

                        if ((mesh1!=mesh2) && mesh2.current_barrier()) // Computing D* block
                            operatorD(mesh1,mesh2,mat_temp,Dcoeff,gauss_order,true);

                        // Computing N block

                        if ((mesh1!=mesh2) || (mesh1!=cortex))
                            operatorN(mesh1,mesh2,mat_temp,Ncoeff,gauss_order);
                    }
                    if (&mesh1==&mesh2)
                        break;
                }

            // Deflate all current barriers as one

            deflate(mat_temp,geo);

            // Copy mat_temp into mat except the lines for cortex vertices [i_vb_c, i_ve_c] and cortex triangles [i_tb_c, i_te_c].

            mat = Matrix(Nl, Nc);
            mat.set(0.0);
            unsigned iNl = 0;
            for (const auto& mesh : geo.meshes())
                if (mesh!=cortex) {
                    for (const auto& vertex : mesh.vertices())
                        mat.setlin(iNl++,mat_temp.getlin(vertex->index()));
                    if (!mesh.current_barrier())
                        for (const auto& triangle : mesh.triangles())
                            mat.setlin(iNl++,mat_temp.getlin(triangle.index()));
                }

            // ** Construct P: the null-space projector **

            Matrix W;
            {
                Matrix U;
                SparseMatrix s;
                mat.svd(U, s, W);
            }

            // Set S to 0 everywhere, except in the last part of the diag:

            SparseMatrix S(Nc,Nc);
            for (unsigned i=Nl;i<Nc;++i)
                S(i,i) = 1.0;

            P = (W.transpose()*S)*W; // P is a projector: P^2 = P and mat*P*X = 0
            if (filename.length()!=0) {
                std::cout << "Saving projector P (" << filename << ")." << std::endl;
                P.save(filename);
            }
        } else {
            std::cout << "Loading projector P (" << filename << ")." << std::endl;
            P.load(filename);
        }

        // ** Get the gradient of P1&P0 elements on the meshes **
        Matrix MM(M.transpose() * M);
        SymMatrix RR(Nc,Nc); RR.set(0.);
        for (const auto& mesh : geo.meshes())
            mesh.gradient_norm2(RR);

        /// Choose Regularization parameter

        SparseMatrix alphas(Nc,Nc); // diagonal matrix
        Matrix Z;
        double alpha1 = alpha;
        double beta1  = beta;
        if (alpha1<0) { // try an automatic method... TODO find better estimation
            const double nRR_v = RR.submat(0,geo.vertices().size(),0,geo.vertices().size()).frobenius_norm();
            alphas.set(0.);
            alpha1 = MM.frobenius_norm()/(1.e3*nRR_v);
            beta1  = alpha1*50000.;
            std::cout << "AUTOMATIC alphas = " << alpha1 << "\tbeta = " << beta1 << std::endl;
        } else {
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

        Matrix rhs = P.transpose()*M.transpose();
        mat = P*Z.pinverse()*rhs;
    }

    CorticalMat2::CorticalMat2(const Geometry& geo,const Head2EEGMat& M,const std::string& domain_name,
                               const unsigned gauss_order,const double gamma,const std::string &filename)
    {
        Matrix& mat = *this;

        // Re-writting of the optimization problem in M. Clerc, J. Kybic "Cortical mapping by Laplace–Cauchy transmission using a boundary element method".
        // with a Lagrangian formulation as in see http://www.math.uh.edu/~rohop/fall_06/Chapter3.pdf eq3.3
        // find argmin(norm(gradient(X)) under constraints: 
        // H * X = 0 and M * X = m
        // let G be the gradient norm matrix, l1, l2 the lagrange parameters
        // 
        // [ G  H' M'] [   X    ]   [ 0 ]
        // | H  0    | |   l1   | = | 0 |
        // [ M     0 ] [   l2   ]   [ m ]
        //
        // {----,----}
        //      K
        // we want a submat of the inverse of K (using blockwise inversion, (TODO maybe iterative solution better ?)).
        // Assumptions:
        // - domain_name: the domain containing the sources is an innermost domain (defined as the interior of only one interface (called Cortex)
        // - Cortex interface is composed of one mesh only (no shared vertices)

        const Domain& SourceDomain = geo.domain(domain_name);
        const Interface& Cortex    = SourceDomain.boundaries().front().interface();
        const Mesh& cortex         = Cortex.begin()->mesh();
        
        om_error(SourceDomain.boundaries().size()==1);
        om_error(Cortex.size()==1);

        // shape of the new matrix:
        unsigned Nl = geo.nb_parameters()-geo.nb_current_barrier_triangles()-Cortex.nb_vertices()-Cortex.nb_triangles();
        unsigned Nc = geo.nb_parameters()-geo.nb_current_barrier_triangles();
        std::fstream f(filename.c_str());
        Matrix H;
        if ( !f ) {
            // build the HeadMat:
            // The following is the same as assemble_HM except N_11, D_11 and S_11 are not computed.
            SymMatrix mat_temp(Nc);
            mat_temp.set(0.0);
            double K = 1.0 / (4.0 * M_PI);
            // We iterate over the meshes (or pair of domains) to fill the lower half of the HeadMat (since its symmetry)

            for (const auto& mesh1 : geo.meshes())
                for (const auto& mesh2 : geo.meshes()) {
                    // if mesh1 and mesh2 communicate, i.e they are used for the definition of a common domain
                    const int orientation = geo.oriented(mesh1, mesh2); // equals  0, if they don't have any domains in common
                    // equals  1, if they are both oriented toward the same domain
                    // equals -1, if they are not
                    if ( orientation != 0) {
                        double Scoeff =   orientation * geo.sigma_inv(mesh1, mesh2) * K;
                        double Dcoeff = - orientation * geo.indicator(mesh1, mesh2) * K;
                        double Ncoeff;
                        if ( !(mesh1.current_barrier() || mesh2.current_barrier()) && ( (mesh1 != mesh2)||( mesh1 != cortex) ) ) {
                            // Computing S block first because it's needed for the corresponding N block
                            operatorS(mesh1, mesh2, mat_temp, Scoeff, gauss_order);
                            Ncoeff = geo.sigma(mesh1, mesh2)/geo.sigma_inv(mesh1, mesh2);
                        } else {
                            Ncoeff = orientation * geo.sigma(mesh1, mesh2) * K;
                        }
                        if ( !mesh1.current_barrier() && (( (mesh1 != mesh2)||( mesh1 != cortex) )) ) {
                            // Computing D block
                            operatorD(mesh1, mesh2, mat_temp, Dcoeff, gauss_order);
                        }
                        if ( ( mesh1 != mesh2 ) && ( !mesh2.current_barrier() ) ) {
                            // Computing D* block
                            operatorD(mesh1, mesh2, mat_temp, Dcoeff, gauss_order, true);
                        }
                        // Computing N block
                        if ( (mesh1 != mesh2)||( mesh1 != cortex) ) {
                            operatorN(mesh1, mesh2, mat_temp, Ncoeff, gauss_order);
                        }
                    }
                    if (&mesh1==&mesh2)
                        break;
                }

            // Deflate all current barriers as one
            deflate(mat_temp,geo);

            H = Matrix(Nl + M.nlin(), Nc);
            H.set(0.0);
            // copy mat_temp into H except the lines for cortex vertices [i_vb_c, i_ve_c] and cortex triangles [i_tb_c, i_te_c].
            unsigned iNl = 0;
            for (const auto& mesh : geo.meshes())
                if (mesh!=cortex) {
                    for (const auto& vertex : mesh.vertices())
                        H.setlin(iNl++,mat_temp.getlin(vertex->index()));
                    if (!mesh.current_barrier())
                        for (const auto& triangle : mesh.triangles())
                            H.setlin(iNl++,mat_temp.getlin(triangle.index()));
                }

            if (filename.length()!=0) {
                std::cout << "Saving matrix H (" << filename << ")." << std::endl;
                H.save(filename);
            }
        } else {
            std::cout << "Loading matrix H (" << filename << ")." << std::endl;
            H.load(filename);
        }

        // Concat M to H

        for (unsigned i=Nl;i<Nl+M.nlin();++i)
            for (unsigned j=0;j<Nc;++j)
                H(i,j) = M(i-Nl,j);

        // Get the gradient of P1&P0 elements on the meshes.

        SymMatrix G(Nc);
        G.set(0.);
        for (const auto& mesh : geo.meshes())
            mesh.gradient_norm2(G);

        // Multiply by gamma the submat of current gradient norm2

        for (const auto& mesh : geo.meshes())
            if (!mesh.current_barrier())
                for (const auto& triangle1 : mesh.triangles())
                    for (const auto& triangle2 : mesh.triangles())
                        G(triangle1.index(),triangle1.index()) *= gamma;

        std::cout << "gamma = " << gamma << std::endl;

        G.invert();
        mat = (G*H.transpose()*(H*G*H.transpose()).inverse()).submat(0,Nc,Nl,M.nlin());
    }

    Surf2VolMat::Surf2VolMat(const Geometry& geo,const Matrix& points) {

        // Find the points per domain and generate the indices for the m_points

        std::map<const Domain*,Vertices> m_points;
        unsigned index = 0;
        for (unsigned i=0;i<points.nlin();++i) {
            const Domain& domain = geo.domain(Vect3(points(i,0),points(i,1),points(i,2))); // TODO: see Vertex below....
            if (domain.conductivity()==0.0) {
                std::cerr << " Surf2Vol: Point [ " << points.getlin(i);
                std::cerr << "] is inside a non-conductive domain. Point is dropped." << std::endl;
            } else {
                m_points[&domain].push_back(Vertex(points(i,0), points(i,1),points(i,2),index++));
            }
        }

        Matrix& mat = *this;
        const double K = 1.0/(4.0*M_PI);

        unsigned size = 0; // total number of inside points
        for (const auto& map_element : m_points)
            size += map_element.second.size();

        mat = Matrix(size,(geo.nb_parameters()-geo.nb_current_barrier_triangles()));
        mat.set(0.0);

        for (const auto& map_element : m_points)
            for (const auto& mesh : geo.meshes()) {
                const int orientation = map_element.first->mesh_orientation(mesh);
                if (orientation!=0) {
                    operatorDinternal(mesh,mat,map_element.second,-orientation*K);
                    if (!mesh.current_barrier())
                        operatorSinternal(mesh,mat,map_element.second,orientation*K/map_element.first->conductivity());
                }
            }
    }
}
