/* FILE: $Id$ */

/*
Project Name : OpenMEEG

version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Maureen.Clerc.AT.sophia.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.sophia.inria.fr)

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

#ifndef OPENMEEG_INVERSE_H
#define OPENMEEG_INVERSE_H

#include "matrix.h"
#include "symmatrix.h"
#include "vector.h"
#include "sparse_matrix.h"
#include "fast_sparse_matrix.h"

#include "DLLDefinesOpenMEEG.h"

#define EPSILON 1e-6
#define MINRES_TOL 1e-5

namespace OpenMEEG {

    // static double k2;
    // struct tv_evaluator;
    //
    // inline double tik (const double &x)
    // {
    //     return 0.5*x*x;
    // }
    // inline double tikp (const double &x)
    // {
    //     return x;
    // }
    //
    // inline double tikpp (const double &x)
    // {
    //     return 1.0;
    // }
    //
    // inline double ftv (const double &x)
    // {
    //     return x;
    // }
    // inline double ftvp (const double &x)
    // {
    //     return 1.0;
    // }
    //
    // inline double ftvpp (const double &x)
    // {
    //     return 0.0;
    // }
    //
    // inline double pm (const double &x)
    // {
    //     return -0.5*k2*(exp(-(x*x)/k2)-1);
    // }
    //
    // inline double pmp (const double &x)
    // {
    //     return x*exp(-(x*x)/k2);
    // }
    //
    // inline double aub (const double &x)
    // {
    //     return (sqrt(1+x*x*k2)-1);
    // }
    //
    // inline double aubp (const double &x)
    // {
    //     return x/(k2*sqrt((k2+x*x)*k2));
    // }
    //
    // inline double aubpp (const double &x)
    // {
    //     return 1.0/((k2+x*x)*sqrt((k2+x*x)*k2));
    // }
    //
    // static double (*ftab[4]) (const double &)={0,tik,pm,aub};
    // static double (*fptab[4]) (const double &)={0,tikp,pmp,aubp};
    // static double (*fpptab[4]) (const double &)={0,tikpp,0,aubpp};

    inline Vector gentv( Vector x,
                          const FastSparseMatrix &mat,
                          const FastSparseMatrix &mat_t,
                          const Vector &Ai, double *tv=NULL,
                          double (*f) (const double &)=0,
                          double (*fp) (const double&)=0 )
    {
        Vector v = mat * x;
        Vector grad_norms( v.size()/3 );
        Vector grad_norms_inv( v.size()/3 );
        for(size_t i=0;i<grad_norms.size();i++) {
            double *pt=&v(3*i);
            grad_norms(i) = sqrt(pt[0]*pt[0]+pt[1]*pt[1]+pt[2]*pt[2]);
            grad_norms_inv(i) = grad_norms(i)!=0 ? 1.0/(grad_norms(i)+EPSILON) : 0;
            double normaliz = grad_norms_inv(i)*Ai(i);
            if (fp!=0) normaliz *= fp(grad_norms(i));
            pt[0]*=normaliz; pt[1]*=normaliz; pt[2]*=normaliz;
        }

        if (tv!=NULL && f!=0) {*tv=0; for(size_t i=0;i<grad_norms.size();i++) *tv+=f(grad_norms(i))*Ai(i);}
        if (tv!=NULL && f==0) {*tv=0; for(size_t i=0;i<grad_norms.size();i++) *tv+=grad_norms(i)*Ai(i);}

        return mat_t*v;
    }

    inline double compute_one_tv(Vector x,
                              const FastSparseMatrix &mat,
                              const Vector &Ai,
                              double (*f) (const double &)=0)
    {
        double tv = 0;
        Vector v = mat * x;
        Vector grad_norms( v.size()/3 );
        for(size_t i=0;i<grad_norms.size();i++) {
            double *pt = &v(3*i);
            grad_norms(i) = sqrt(pt[0]*pt[0]+pt[1]*pt[1]+pt[2]*pt[2]);
            if (f!=0) {
                tv += f(grad_norms(i))*Ai(i);
            } else {
                tv += grad_norms(i)*Ai(i);
            }
        }
        return tv;
    }

    // ========================================================
    // = Define Hessian matrices for linear inversion methods =
    // ========================================================

    class OPENMEEG_EXPORT MN_Hessian {

        const Matrix &Transfer;
        const double alpha;

    public:

        MN_Hessian(const Matrix &TransferMat, const double &Alpha):Transfer(TransferMat),alpha(Alpha) {}

        inline Vector operator * (const Vector &x) const { return Transfer.tmult(Transfer*x)+alpha*x; }
    };

    class OPENMEEG_EXPORT WMN_Hessian {

        const Matrix &Transfer;
        const double alpha;
        Vector weights;

    public:

        WMN_Hessian(const Matrix &TransferMat, const double &Alpha):Transfer(TransferMat),alpha(Alpha) {
            Vector v(Transfer.ncol());
            for(size_t i = 0; i < weights.size(); ++i) {
                Vector col = Transfer.getcol(i);
                v(i) = pow(col.norm(),2);
            }
            weights = v;
        }

        inline Vector operator * (const Vector &x) const {
            return Transfer.tmult(Transfer*x)+alpha*(weights.kmult(x));
        }

    };

    class OPENMEEG_EXPORT HEAT_Hessian {

        const Matrix &m_transfer;
        const FastSparseMatrix &m_mat;
        const FastSparseMatrix &m_mat_t;
        const double m_alpha;

    public:

        HEAT_Hessian(const Matrix &transfer,
                         const FastSparseMatrix &mat,
                         const FastSparseMatrix &mat_t,
                         const double &alpha):
        m_transfer(transfer),m_mat(mat),m_mat_t(mat_t),m_alpha(alpha) {}

        inline Vector operator * ( const Vector &x) const {
            return m_transfer.tmult(m_transfer*x)+m_alpha*(m_mat_t*(m_mat*x));
        }
    };

    // ========================================================

    template<class T> // T should be a linear operator
    size_t MinRes2(const T& A,const Vector& b,Vector& x0,double tol) {

        size_t n_max=10000;
        size_t n=1; size_t N=x0.size();
        Vector v(N); v.set(0.0);
        Vector v_hat=b-A*x0;
        double beta=v_hat.norm();
        Vector v_old(v.size());
        Vector Av(v.size());
        double c=1; double c_old=1; double s_old=0; double s=0;
        Vector w(N); w.set(0.0);
        Vector w_oold(N); Vector w_old(w,DEEP_COPY);
        double eta=beta;
        Vector xMR=x0;
        double norm_rMR=beta; double norm_r0=beta;
        double c_oold,s_oold,r1_hat,r1,r2,r3,alpha,beta_old;
        while ((n < n_max+1) && (norm_rMR/norm_r0 > tol) ) {
            n=n+1;
            //Lanczos
            v_old=v;
            v=v_hat*(1.0/beta); Av=A*v; alpha=v*Av;
            v_hat=Av-alpha*v-beta*v_old;
            beta_old=beta; beta=v_hat.norm();
            //QR factorization
            c_oold=c_old; c_old=c;  s_oold=s_old; s_old=s;
            r1_hat=c_old*alpha-c_oold*s_old*beta_old;
            r1 = sqrt(r1_hat*r1_hat+beta*beta);
            r2 = s_old*alpha+c_oold*c_old*beta_old;
            r3 = s_oold*beta_old;
            //Givens rotation
            c=r1_hat/r1;
            s=beta/r1;
            //update
            w_oold=w_old; w_old=w;
            w=(v-r3*w_oold-r2*w_old)*(1.0/r1);
            xMR+=c*eta*w; norm_rMR=norm_rMR*fabs(s);
            eta=-s*eta;
        }
        std::cout<<"\r";
        return n;
    }

    // ===========================================
    // = Define all the linear inversion methods =
    // ===========================================

    template<class T>
    void LIN_inverse (Matrix& EstimatedData, const T& hess, const Matrix& GainMatrix, const Matrix& Data) {
        size_t nT = Data.ncol();
        EstimatedData = Matrix(GainMatrix.ncol(),nT);

        #ifdef USE_OMP
        #pragma omp parallel for
        #endif
        for(long frame=0;frame<(long)nT;frame++) { // loop over frame
            Vector m = Data.getcol((size_t) frame);
            Vector v(GainMatrix.ncol()); v.set(0.0);

            //==========  Invert =======================//
            size_t niter = MinRes2(hess,GainMatrix.tmult(m),v,MINRES_TOL);

            for(size_t i=0;i<EstimatedData.nlin();i++) EstimatedData(i,(size_t) frame) = v(i);

            #ifdef USE_OMP
            #pragma omp critical
            #endif
            std::cout << ">> Frame " << frame+1 << " / " << nT
                      << " : Rel. Err. = " << (GainMatrix*v-m).norm()/m.norm()
                      << " : Nb. iter. MinRes = " << niter
                      << std::endl;
        }
    }

    void compute_mn (Matrix& EstimatedData, const Matrix& Data, const Matrix& GainMatrix, double SmoothWeight) {
        Matrix eye(GainMatrix.nlin(),GainMatrix.nlin());
        eye.set(0);
        for (size_t i=0;i<GainMatrix.nlin();++i)
            eye(i,i) = SmoothWeight;
        EstimatedData = GainMatrix.transpose() * (GainMatrix * GainMatrix.transpose() + eye).inverse() * Data;
    }

    // ================= Iterative Mininum norm inversion =======================//

    class IMN_inverse: public Matrix {
    public:
        IMN_inverse (const Matrix& Data, const Matrix& GainMatrix, double SmoothWeight);
        virtual ~IMN_inverse () {};
    };

    IMN_inverse::IMN_inverse (const Matrix& Data, const Matrix& GainMatrix, double SmoothWeight) {
        std::cout << "Running Iterative MN inversion" << std::endl;
        MN_Hessian hess(GainMatrix,SmoothWeight);
        LIN_inverse(*this,hess,GainMatrix,Data);
    }

    // ================= Mininum norm inversion =======================//

    class MN_inverse: public Matrix {
    public:
        MN_inverse (const Matrix& Data, const Matrix& GainMatrix, double SmoothWeight);
        virtual ~MN_inverse () {};
    };

    MN_inverse::MN_inverse (const Matrix& Data, const Matrix& GainMatrix, double SmoothWeight) {
        std::cout << "Running MN inversion" << std::endl;
        compute_mn(*this,Data,GainMatrix,SmoothWeight);
    }

    // ================= Weighted Mininum norm inversion =======================//

    class WMN_inverse: public Matrix {
    public:
        WMN_inverse (const Matrix& Data, const Matrix& GainMatrix, double SmoothWeight);
        virtual ~WMN_inverse () {};
    };

    WMN_inverse::WMN_inverse (const Matrix& Data, const Matrix& GainMatrix, double SmoothWeight) {
        std::cout << "Running WMN inversion" << std::endl;
        WMN_Hessian hess(GainMatrix,SmoothWeight);
        LIN_inverse(*this,hess,GainMatrix,Data);
    }

    // ================= Gradient based Mininum norm inversion ================ //

    class HEAT_inverse: public Matrix {
    public:
        HEAT_inverse (const Matrix& Data, const Matrix& GainMatrix, const SparseMatrix& SmoothMatrix, double SmoothWeight);
        virtual ~HEAT_inverse () {};
    };

    HEAT_inverse::HEAT_inverse (const Matrix& Data, const Matrix& GainMatrix, const SparseMatrix& SmoothMatrix, double SmoothWeight) {
        std::cout << "Running HEAT inversion" << std::endl;
        FastSparseMatrix fastSmoothMatrix(SmoothMatrix);
        FastSparseMatrix fastSmoothMatrix_t(SmoothMatrix.transpose());
        HEAT_Hessian hess(GainMatrix,fastSmoothMatrix,fastSmoothMatrix_t,SmoothWeight);
        LIN_inverse(*this,hess,GainMatrix,Data);
    }

    // ================= Total variation based inversion =================== //

    void compute_tv(Matrix& EstimatedData, const Matrix& Data, const Matrix& GainMatrix, const SparseMatrix& SmoothMatrix, const Vector& AiVector, double SmoothWeight, size_t MaxNbIter, double StoppingTol)
    {
        FastSparseMatrix fastSmoothMatrix(SmoothMatrix);
        FastSparseMatrix fastSmoothMatrix_t(SmoothMatrix.transpose());

        size_t nT = Data.ncol();
        EstimatedData = Matrix(GainMatrix.ncol(),nT);

        // #ifdef USE_OMP
        // #pragma omp parallel for
        // #endif
        for(size_t frame=0;frame<nT;frame++) {
            std::cout << ">> Frame " << frame+1 << " / " << nT << std::endl;
            Vector m = Data.getcol(frame);
            Vector v(EstimatedData.nlin());

            // ====================  initialization of source vector ===================== //
            if(frame==0) v.set(0.0);
            else v = EstimatedData.getcol(frame-1);

            double tv_v = compute_one_tv(v,fastSmoothMatrix,AiVector);

            bool errorTest = true;

            // ========  Backtracking line search parameters for gradient step  ========= //
            double alpha = 0.001;
            double beta = 0.5;
            int max_iter_line_search = 10;

            // ================== Inverse problem via gradient descent ================== //
            size_t t;
            for(t=0;t<MaxNbIter && errorTest;t++) {
                Vector gradtv = gentv(v,fastSmoothMatrix,fastSmoothMatrix_t,AiVector);
                Vector err_vec = GainMatrix*v-m;
                Vector graddata = GainMatrix.tmult(err_vec);
                Vector grad = (-SmoothWeight)*gradtv - graddata;
                double f_v = pow(err_vec.norm(),2) + SmoothWeight*tv_v;

                // ======= Backtracking line search for gradient step ======= //
                double search_slope = alpha*grad.norm();
                double f_v_dv;
                double tv_v_dv;
                Vector v_dv;

                double grad_step = 1.0;

    #define USE_LINE_SEARCH 1

    #if USE_LINE_SEARCH
                int iter_line_search = 0;
                bool stop_line_search = false;
                while ( stop_line_search != true && (++iter_line_search < max_iter_line_search) ) {
                    v_dv = v+grad_step*grad;
                    double f_v_dv_data = pow((m-GainMatrix*(v_dv)).norm(),2);
                    tv_v_dv = compute_one_tv(v_dv,fastSmoothMatrix,AiVector);
                    f_v_dv = f_v_dv_data + SmoothWeight*tv_v_dv;
                    if ( grad_step*search_slope < (f_v - f_v_dv)) {
                        stop_line_search = true;
                    } else {
                        grad_step = beta*grad_step;
                    }
                }
    #elif
                v_dv = v+grad_step*grad;
    #endif

                double tol = (v_dv-v).norm()/v.norm();

                v = v_dv;
                tv_v = tv_v_dv;
                errorTest = tol>StoppingTol;// || iter_line_search<max_iter_line_search;

                if ((t%100)==0 || !errorTest || (t == (MaxNbIter-1)))
                    printf("Energy %e   Relative Error %f   TV %f   Tol %e   GradStep %f Iter %d\n",
                           f_v,(err_vec).norm()/m.norm(),tv_v,tol,grad_step,(int)t);
            }
            printf("Total number of iterations : %d\n",(int)t);

            //===========================================================================//
            EstimatedData.setcol(frame,v);
        }
    }

    class TV_inverse: public Matrix {
    public:
        TV_inverse (const Matrix& Data, const Matrix& GainMatrix, const SparseMatrix& SmoothMatrix, const Vector& AiVector, double SmoothWeight, size_t MaxNbIter, double StoppingTol);
        virtual ~TV_inverse () {};
    };

    TV_inverse::TV_inverse (const Matrix& Data, const Matrix& GainMatrix, const SparseMatrix& SmoothMatrix, const Vector& AiVector, double SmoothWeight, size_t MaxNbIter, double StoppingTol) {
        std::cout << "Running TV inversion" << std::endl;
        compute_tv(*this,Data,GainMatrix,SmoothMatrix,AiVector,SmoothWeight,MaxNbIter,StoppingTol);
    }
}

#endif /* OPENMEEG_INVERSE_H */

