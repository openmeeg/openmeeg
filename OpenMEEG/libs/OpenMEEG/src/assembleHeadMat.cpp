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

    template<class T>
    void deflat(T& M, const Interface& i, double coef) 
    {
        // deflate the Matrix
        for (Interface::const_iterator omit=i.begin();omit!=i.end();++omit) {
            for (Mesh::const_vertex_iterator vit1=omit->mesh().vertex_begin();vit1!=omit->mesh().vertex_end();++vit1) {
                #pragma omp parallel for
                #ifndef OPENMP_3_0
                for (int i2=vit1-omit->mesh().vertex_begin();i2<omit->mesh().vertex_size();++i2) {
                    const Mesh::const_vertex_iterator vit2 = omit->mesh().vertex_begin()+i2;
                #else
                for (Mesh::const_vertex_iterator vit2=vit1;vit2<omit->mesh().vertex_end();++vit2) {
                #endif
                    M((*vit1)->index(),(*vit2)->index()) += coef;
                }
            }
        }
    }

    template<class T>
    void deflat(T& M, const Geometry& geo)
    {
        //deflat all current barriers as one
        unsigned nb_vertices=0,i_first=0;
        double coef=0.0;
        for ( std::vector<Strings>::const_iterator git = geo.geo_group().begin(); git != geo.geo_group().end(); ++git) {
            nb_vertices=0;
            i_first=0;
            for(std::vector<std::string>::const_iterator mit=git->begin();mit!=git->end();++mit){
                const Mesh msh=geo.mesh(*mit);
                if(msh.outermost()){
                    nb_vertices+=msh.nb_vertices();
                    if(i_first==0)
                        i_first=(*msh.vertex_begin())->index();
                }
            }
            coef=M(i_first,i_first)/nb_vertices;
            for(std::vector<std::string>::const_iterator mit=git->begin();mit!=git->end();++mit){
                Mesh msh=geo.mesh(*mit);
                if(msh.outermost())
                    for(Mesh::const_vertex_iterator vit1=msh.vertex_begin();vit1!=msh.vertex_end();++vit1){
                        #pragma omp parallel for
                        #ifndef OPENMP_3_0
                        for (int i2=vit1-msh.vertex_begin();i2<msh.vertex_size();++i2) {
                            const Mesh::const_vertex_iterator vit2 = msh.vertex_begin()+i2;
                        #else
                        for (Mesh::const_vertex_iterator vit2=vit1;vit2<msh.vertex_end();++vit2) {
                        #endif
                            M((*vit1)->index(),(*vit2)->index()) += coef;
                        }
                    }
            }
        }
    }

    void assemble_HM(const Geometry& geo, SymMatrix& mat, const unsigned gauss_order) 
    {
        mat = SymMatrix((geo.size()-geo.nb_current_barrier_triangles()));
        mat.set(0.0);
        double K = 1.0 / (4.0 * M_PI);

        // We iterate over the meshes (or pair of domains) to fill the lower half of the HeadMat (since its symmetry)
        for(Geometry::const_iterator mit1 = geo.begin(); mit1 != geo.end(); ++mit1) {
            if(!mit1->isolated()){
                for(Geometry::const_iterator mit2 = geo.begin(); (mit2 != (mit1+1)); ++mit2) {
                    if((!mit2->isolated()) && (geo.sigma(*mit1,*mit2)!=0.0)){
                        // if mit1 and mit2 communicate, i.e they are used for the definition of a common domain
                        const int orientation = geo.oriented(*mit1, *mit2); // equals  0, if they don't have any domains in common
                                                                            // equals  1, if they are both oriented toward the same domain
                                                                            // equals -1, if they are not
                        if(orientation!=0){
                            double Scoeff =   orientation * geo.sigma_inv(*mit1, *mit2) * K;
                            double Dcoeff = - orientation * geo.indicator(*mit1, *mit2) * K;
                            double Ncoeff;

                            if( (!mit1->current_barrier()) && (!mit2->current_barrier()) ) {
                                // Computing S block first because it's needed for the corresponding N block
                                operatorS(*mit1, *mit2, mat, Scoeff, gauss_order);
                                Ncoeff = geo.sigma(*mit1, *mit2)/geo.sigma_inv(*mit1, *mit2);
                            }else{
                                Ncoeff = orientation * geo.sigma(*mit1, *mit2) * K;
                            }

                            if(!mit1->current_barrier()){
                                // Computing D block
                                operatorD(*mit1, *mit2, mat, Dcoeff, gauss_order,false);
                            }
                            if((*mit1!=*mit2) && (!mit2->current_barrier())){
                                // Computing D* block
                                operatorD(*mit1, *mit2, mat, Dcoeff, gauss_order, true);
                            }
                            // Computing N block
                            operatorN(*mit1, *mit2, mat, Ncoeff, gauss_order);
                        }
                    }
                }
            }
        }
        // Deflate all current barriers as one
        deflat(mat,geo);
    }

    void assemble_cortical(const Geometry& geo, Matrix& mat, const Head2EEGMat& M, const std::string& domain_name, const unsigned gauss_order, double alpha, double beta, const std::string &filename)
    {
        // Following the article: M. Clerc, J. Kybic "Cortical mapping by Laplace–Cauchy transmission using a boundary element method".
        // Assumptions:
        // - domain_name: the domain containing the sources is an innermost domain (defined as the interior of only one interface (called Cortex))
        // - Cortex interface is composed of one mesh only (no shared vertices)

        const Domain& SourceDomain  = geo.domain(domain_name);
        const Interface& Cortex     = SourceDomain.begin()->interface();
        const Mesh& cortex          = Cortex.begin()->mesh();
        
        om_error(SourceDomain.size()==1);
        om_error(Cortex.size()==1);

        // shape of the new matrix:
        unsigned Nl = geo.size()-geo.nb_current_barrier_triangles()-Cortex.nb_vertices()-Cortex.nb_triangles();
        unsigned Nc = geo.size()-geo.nb_current_barrier_triangles();
        std::fstream f(filename.c_str());
        Matrix P;
        if ( !f ) {
            // build the HeadMat:
            // The following is the same as assemble_HM except N_11, D_11 and S_11 are not computed.
            SymMatrix mat_temp(Nc);
            mat_temp.set(0.0);
            double K = 1.0 / (4.0 * M_PI);
            // We iterate over the meshes (or pair of domains) to fill the lower half of the HeadMat (since its symmetry)
            for ( Geometry::const_iterator mit1 = geo.begin(); mit1 != geo.end(); ++mit1) {
                for ( Geometry::const_iterator mit2 = geo.begin(); (mit2 != (mit1+1)); ++mit2) {
                    // if mit1 and mit2 communicate, i.e they are used for the definition of a common domain
                    const int orientation = geo.oriented(*mit1, *mit2); // equals  0, if they don't have any domains in common
                    // equals  1, if they are both oriented toward the same domain
                    // equals -1, if they are not
                    if ( orientation != 0) {
                        double Scoeff =   orientation * geo.sigma_inv(*mit1, *mit2) * K;
                        double Dcoeff = - orientation * geo.indicator(*mit1, *mit2) * K;
                        double Ncoeff;
                        if ( !(mit1->current_barrier() || mit2->current_barrier()) && ( (*mit1 != *mit2)||( *mit1 != cortex) ) ) {
                            // Computing S block first because it's needed for the corresponding N block
                            operatorS(*mit1, *mit2, mat_temp, Scoeff, gauss_order);
                            Ncoeff = geo.sigma(*mit1, *mit2)/geo.sigma_inv(*mit1, *mit2);
                        } else {
                            Ncoeff = orientation * geo.sigma(*mit1, *mit2) * K;
                        }
                        if ( !mit1->current_barrier() && (( (*mit1 != *mit2)||( *mit1 != cortex) )) ) {
                            // Computing D block
                            operatorD(*mit1, *mit2, mat_temp, Dcoeff, gauss_order,false);
                        }
                        if ( ( *mit1 != *mit2 ) && ( !mit2->current_barrier() ) ) {
                            // Computing D* block
                            operatorD(*mit1, *mit2, mat_temp, Dcoeff, gauss_order, true);
                        }
                        // Computing N block
                        if ( (*mit1 != *mit2)||( *mit1 != cortex) ) {
                            operatorN(*mit1, *mit2, mat_temp, Ncoeff, gauss_order);
                        }
                    }
                }
            }
            // Deflate all current barriers as one
            deflat(mat_temp,geo);

            mat = Matrix(Nl, Nc);
            mat.set(0.0);
            // copy mat_temp into mat except the lines for cortex vertices [i_vb_c, i_ve_c] and cortex triangles [i_tb_c, i_te_c].
            unsigned iNl = 0;
            for ( Geometry::const_iterator mit = geo.begin(); mit != geo.end(); ++mit) {
                if ( *mit != cortex ) {
                    for ( Mesh::const_vertex_iterator vit = mit->vertex_begin(); vit != mit->vertex_end(); ++vit) {
                        mat.setlin(iNl, mat_temp.getlin((*vit)->index()));
                        ++iNl;
                    }
                    if ( !mit->current_barrier() ) {
                        for ( Mesh::const_iterator tit = mit->begin(); tit != mit->end(); ++tit) {
                            mat.setlin(iNl, mat_temp.getlin(tit->index()));
                            ++iNl;
                        }
                    }
                }
            }
            // ** Construct P: the null-space projector **
            Matrix W;
            {
                Matrix U;
                SparseMatrix s;
                mat.svd(U, s, W);
            }

            SparseMatrix S(Nc,Nc);
            // we set S to 0 everywhere, except in the last part of the diag:
            for ( unsigned i = Nl; i < Nc; ++i) {
                S(i, i) = 1.0;
            }
            P = (W.transpose() * S) * W; // P is a projector: P^2 = P and mat*P*X = 0
            if ( filename.length() != 0 ) {
                std::cout << "Saving projector P (" << filename << ")." << std::endl;
                P.save(filename);
            }
        } else {
            std::cout << "Loading projector P (" << filename << ")." << std::endl;
            P.load(filename);
        }

        // ** Get the gradient of P1&P0 elements on the meshes **
        Matrix MM(M.transpose() * M);
        SymMatrix RR(Nc, Nc); RR.set(0.);
        for ( Geometry::const_iterator mit = geo.begin(); mit != geo.end(); ++mit) {
            mit->gradient_norm2(RR);
        }

        // ** Choose Regularization parameter **
        SparseMatrix alphas(Nc,Nc); // diagonal matrix
        Matrix Z;
        if ( alpha < 0 ) { // try an automatic method... TODO find better estimation
            double nRR_v = RR.submat(0, geo.nb_vertices(), 0, geo.nb_vertices()).frobenius_norm();
            alphas.set(0.);
            alpha = MM.frobenius_norm() / (1.e3*nRR_v);
            beta  = alpha * 50000.;
            for ( Vertices::const_iterator vit = geo.vertex_begin(); vit != geo.vertex_end(); ++vit) {
                alphas(vit->index(), vit->index()) = alpha;
            }
            for ( Meshes::const_iterator mit = geo.begin(); mit != geo.end(); ++mit) {
                if ( !mit->current_barrier() ) {
                    for ( Mesh::const_iterator tit = mit->begin(); tit != mit->end(); ++tit) {
                        alphas(tit->index(), tit->index()) = beta;
                    }
                }
            }
            std::cout << "AUTOMATIC alphas = " << alpha << "\tbeta = " << beta << std::endl;
        } else {
            for ( Vertices::const_iterator vit = geo.vertex_begin(); vit != geo.vertex_end(); ++vit) {
                alphas(vit->index(), vit->index()) = alpha;
            }
            for ( Meshes::const_iterator mit = geo.begin(); mit != geo.end(); ++mit) {
                if ( !mit->current_barrier() ) {
                    for ( Mesh::const_iterator tit = mit->begin(); tit != mit->end(); ++tit) {
                        alphas(tit->index(), tit->index()) = beta;
                    }
                }
            }
            std::cout << "alphas = " << alpha << "\tbeta = " << beta << std::endl;
        }
        Z = P.transpose() * (MM + alphas*RR) * P;

        // ** PseudoInverse and return **
        // X = P * { (M*P)' * (M*P) + (R*P)' * (R*P) }¡(-1) * (M*P)'m
        // X = P * { P'*M'*M*P + P'*R'*R*P }¡(-1) * P'*M'm
        // X = P * { P'*(MM + a*RR)*P }¡(-1) * P'*M'm
        // X = P * Z¡(-1) * P' * M'm
        Matrix rhs = P.transpose() * M.transpose();
        mat = P * Z.pinverse() * rhs;
    }

    void assemble_cortical2(const Geometry& geo, Matrix& mat, const Head2EEGMat& M, const std::string& domain_name, const unsigned gauss_order, double gamma, const std::string &filename)
    {
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
        const Interface& Cortex    = SourceDomain.begin()->interface();
        const Mesh& cortex         = Cortex.begin()->mesh();
        
        om_error(SourceDomain.size()==1);
        om_error(Cortex.size()==1);

        // shape of the new matrix:
        unsigned Nl = geo.size()-geo.nb_current_barrier_triangles()-Cortex.nb_vertices()-Cortex.nb_triangles();
        unsigned Nc = geo.size()-geo.nb_current_barrier_triangles();
        std::fstream f(filename.c_str());
        Matrix H;
        if ( !f ) {
            // build the HeadMat:
            // The following is the same as assemble_HM except N_11, D_11 and S_11 are not computed.
            SymMatrix mat_temp(Nc);
            mat_temp.set(0.0);
            double K = 1.0 / (4.0 * M_PI);
            // We iterate over the meshes (or pair of domains) to fill the lower half of the HeadMat (since its symmetry)
            for ( Geometry::const_iterator mit1 = geo.begin(); mit1 != geo.end(); ++mit1) {
                for ( Geometry::const_iterator mit2 = geo.begin(); (mit2 != (mit1+1)); ++mit2) {
                    // if mit1 and mit2 communicate, i.e they are used for the definition of a common domain
                    const int orientation = geo.oriented(*mit1, *mit2); // equals  0, if they don't have any domains in common
                    // equals  1, if they are both oriented toward the same domain
                    // equals -1, if they are not
                    if ( orientation != 0) {
                        double Scoeff =   orientation * geo.sigma_inv(*mit1, *mit2) * K;
                        double Dcoeff = - orientation * geo.indicator(*mit1, *mit2) * K;
                        double Ncoeff;
                        if ( !(mit1->current_barrier() || mit2->current_barrier()) && ( (*mit1 != *mit2)||( *mit1 != cortex) ) ) {
                            // Computing S block first because it's needed for the corresponding N block
                            operatorS(*mit1, *mit2, mat_temp, Scoeff, gauss_order);
                            Ncoeff = geo.sigma(*mit1, *mit2)/geo.sigma_inv(*mit1, *mit2);
                        } else {
                            Ncoeff = orientation * geo.sigma(*mit1, *mit2) * K;
                        }
                        if ( !mit1->current_barrier() && (( (*mit1 != *mit2)||( *mit1 != cortex) )) ) {
                            // Computing D block
                            operatorD(*mit1, *mit2, mat_temp, Dcoeff, gauss_order);
                        }
                        if ( ( *mit1 != *mit2 ) && ( !mit2->current_barrier() ) ) {
                            // Computing D* block
                            operatorD(*mit1, *mit2, mat_temp, Dcoeff, gauss_order, true);
                        }
                        // Computing N block
                        if ( (*mit1 != *mit2)||( *mit1 != cortex) ) {
                            operatorN(*mit1, *mit2, mat_temp, Ncoeff, gauss_order);
                        }
                    }
                }
            }
            // Deflate all current barriers as one
            deflat(mat_temp,geo);

            H = Matrix(Nl + M.nlin(), Nc);
            H.set(0.0);
            // copy mat_temp into H except the lines for cortex vertices [i_vb_c, i_ve_c] and cortex triangles [i_tb_c, i_te_c].
            unsigned iNl = 0;
            for ( Geometry::const_iterator mit = geo.begin(); mit != geo.end(); ++mit) {
                if ( *mit != cortex ) {
                    for ( Mesh::const_vertex_iterator vit = mit->vertex_begin(); vit != mit->vertex_end(); ++vit) {
                        H.setlin(iNl, mat_temp.getlin((*vit)->index()));
                        ++iNl;
                    }
                    if ( !mit->current_barrier() ) {
                        for ( Mesh::const_iterator tit = mit->begin(); tit != mit->end(); ++tit) {
                            H.setlin(iNl, mat_temp.getlin(tit->index()));
                            ++iNl;
                        }
                    }
                }
            }
            if ( filename.length() != 0 ) {
                std::cout << "Saving matrix H (" << filename << ")." << std::endl;
                H.save(filename);
            }
        } else {
            std::cout << "Loading matrix H (" << filename << ")." << std::endl;
            H.load(filename);
        }

        // concat M to H
        for ( unsigned i = Nl; i < Nl + M.nlin(); ++i) {
            for ( unsigned j = 0; j < Nc; ++j) {
                H(i, j) = M(i-Nl, j);
            }
        }

        // ** Get the gradient of P1&P0 elements on the meshes **
        SymMatrix G(Nc);
        G.set(0.);
        for ( Geometry::const_iterator mit = geo.begin(); mit != geo.end(); ++mit) {
            mit->gradient_norm2(G);
        }
        // multiply by gamma the submat of current gradient norm2
        for ( Meshes::const_iterator mit = geo.begin(); mit != geo.end(); ++mit) {
            if ( !mit->current_barrier() ) {
                for ( Mesh::const_iterator tit1 = mit->begin(); tit1 != mit->end(); ++tit1) {
                    for ( Mesh::const_iterator tit2 = mit->begin(); tit2 != mit->end(); ++tit2) {
                        G(tit1->index(), tit2->index()) *= gamma;
                    }
                }
            }
        }
        std::cout << "gamma = " << gamma << std::endl;

        G.invert();
        mat = (G * H.transpose() * (H * G * H.transpose()).inverse()).submat(0, Nc, Nl, M.nlin());
    }

    void assemble_Surf2Vol(const Geometry& geo, Matrix& mat, const std::map<const Domain, Vertices> m_points) 
    {
        const double K = 1.0/(4.0*M_PI);

        unsigned size = 0; // total number of inside points
        for ( std::map<const Domain, Vertices>::const_iterator dvit = m_points.begin(); dvit != m_points.end(); ++dvit) {
            size += dvit->second.size();
        }

        mat = Matrix(size, (geo.size() - geo.nb_current_barrier_triangles()));
        mat.set(0.0);

        for ( std::map<const Domain, Vertices>::const_iterator dvit = m_points.begin(); dvit != m_points.end(); ++dvit) {
            for ( Geometry::const_iterator mit = geo.begin(); mit != geo.end(); ++mit) {
                int orientation = dvit->first.mesh_orientation(*mit);
                if ( orientation != 0 ) {
                    operatorDinternal(*mit, mat, dvit->second, orientation * -1. * K);
                    if ( !mit->current_barrier() ) {
                        operatorSinternal(*mit, mat, dvit->second, orientation * K / geo.sigma(dvit->first));
                    }
                }
            }
        }
    }

    HeadMat::HeadMat(const Geometry& geo, const unsigned gauss_order)
    {
        assemble_HM(geo, *this, gauss_order);
    }

    CorticalMat::CorticalMat(const Geometry& geo, const Head2EEGMat& M, const std::string& domain_name, const unsigned gauss_order, double a, double b, const std::string &filename)
    {
        assemble_cortical(geo, *this, M, domain_name, gauss_order, a, b, filename);
    }

    CorticalMat2::CorticalMat2(const Geometry& geo, const Head2EEGMat& M, const std::string& domain_name, const unsigned gauss_order, double gamma, const std::string &filename)
    {
        assemble_cortical2(geo, *this, M, domain_name, gauss_order, gamma, filename);
    }

    Surf2VolMat::Surf2VolMat(const Geometry& geo, const Matrix& points) 
    {
        std::map<const Domain, Vertices> m_points;

        unsigned index = 0;
        // Find the points per domain and generate the indices for the m_points
        for ( unsigned i = 0; i < points.nlin(); ++i) {
            const Domain domain = geo.domain(Vect3(points(i, 0), points(i, 1), points(i, 2)));
            if ( domain.sigma()==0.0 ) {
                std::cerr << " Surf2Vol: Point [ " << points.getlin(i);
                std::cerr << "] is inside a nonconductive domain. Point is dropped." << std::endl;
            } else {
                m_points[domain].push_back(Vertex(points(i, 0), points(i, 1), points(i, 2), index++));
            }
        }

        assemble_Surf2Vol(geo, *this, m_points);
    }
} // namespace OpenMEEG
