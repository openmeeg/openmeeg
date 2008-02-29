/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
*/

#include "matrice.h"
#include "symmatrice.h"
#include "vecteur.h"
#include "sparse_matrice.h"
#include "fast_sparse_matrice.h"

#define EPSILON 1e-6
#define MINRES_TOL 1e-5

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

inline vecteur gentv( vecteur x,
                      const fast_sparse_matrice &mat,
                      const fast_sparse_matrice &mat_t,
                      const vecteur &Ai, double *tv=NULL,
                      double (*f) (const double &)=0,
                      double (*fp) (const double&)=0 )
{
    vecteur v = mat * x;
    vecteur grad_norms( v.size()/3 );
    vecteur grad_norms_inv( v.size()/3 );
    for(size_t i=0;i<grad_norms.size();i++)
    {
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

inline double compute_tv(vecteur x,
                          const fast_sparse_matrice &mat,
                          const fast_sparse_matrice &mat_t,
                          const vecteur &Ai,
                          double (*f) (const double &)=0,
                          double (*fp) (const double&)=0)
{
    double tv = 0;
    vecteur v = mat * x;
    vecteur grad_norms( v.size()/3 );
    for(size_t i=0;i<grad_norms.size();i++)
    {
        double *pt=&v(3*i);
        grad_norms(i)=sqrt(pt[0]*pt[0]+pt[1]*pt[1]+pt[2]*pt[2]);
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
class MN_Hessian : public LinOp
{
    const matrice &Transfer;
    const double alpha;

public:

    MN_Hessian(const matrice &TransferMat, const double &Alpha):Transfer(TransferMat),alpha(Alpha) {}

    virtual vecteur operator * (const vecteur &x) const
    {
        return Transfer.tmult(Transfer*x)+alpha*x;
    }

};

class WMN_Hessian : public LinOp
{
    const matrice &Transfer;
    const double alpha;
    vecteur weights;

public:

    WMN_Hessian(const matrice &TransferMat, const double &Alpha):Transfer(TransferMat),alpha(Alpha) {
        vecteur v(Transfer.ncol());
        for(size_t i = 0; i < weights.size(); ++i)
        {
            vecteur col = Transfer.getcol(i);
            v(i) = pow(col.norm(),2);
        }
        weights = v;
    }

    virtual vecteur operator * (const vecteur &x) const
    {
        return Transfer.tmult(Transfer*x)+alpha*(weights.kmult(x));
    }

};

class HEAT_Hessian : public LinOp
{
    const matrice &m_transfer;
    const fast_sparse_matrice &m_mat;
    const fast_sparse_matrice &m_mat_t;
    const double m_alpha;
public:
    HEAT_Hessian(const matrice &transfer,
                     const fast_sparse_matrice &mat,
                     const fast_sparse_matrice &mat_t,
                     const double &alpha):
    m_transfer(transfer),m_mat(mat),m_mat_t(mat_t),m_alpha(alpha) {}
    virtual vecteur operator * ( const vecteur &x) const
    {
        return m_transfer.tmult(m_transfer*x)+m_alpha*(m_mat_t*(m_mat*x));
    }
};

// ========================================================

size_t MinRes2(const LinOp& A,const vecteur& b,vecteur& x0,double tol)
{
    size_t n_max=10000;
    size_t n=1; size_t N=x0.size();
    vecteur v(N); v.set(0.0);
    vecteur v_hat=b-A*x0;
    double beta=v_hat.norm();
    vecteur v_old(v.size());
    vecteur Av(v.size());
    double c=1; double c_old=1; double s_old=0; double s=0;
    vecteur w(N); w.set(0.0);
    vecteur w_oold(N); vecteur w_old=w.duplicate();
    double eta=beta;
    vecteur xMR=x0;
    double norm_rMR=beta; double norm_r0=beta;
    double c_oold,s_oold,r1_hat,r1,r2,r3,alpha,beta_old;
    while ((n < n_max+1) && (norm_rMR/norm_r0 > tol) )
    {
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
    }// end while
    std::cout<<"\r";
    return n;
}

// ===========================================
// = Define all the linear inversion methods =
// ===========================================

void LIN_inverse (matrice& EstimatedData, const LinOp& hess, const matrice& GainMatrix, const matrice& Data) {
    size_t nT = Data.ncol();
    EstimatedData = matrice(GainMatrix.ncol(),nT);

    #ifdef USE_OMP
    #pragma omp parallel for
    #endif
    for(int frame=0;frame<nT;frame++)// loop over frame
    {
        vecteur m = Data.getcol(frame);
        vecteur v(GainMatrix.ncol()); v.set(0.0);

        //==========  Invert =======================//
        size_t niter = MinRes2(hess,GainMatrix.tmult(m),v,MINRES_TOL);

        for(size_t i=0;i<EstimatedData.nlin();i++) EstimatedData(i,frame) = v(i);

        #ifdef USE_OMP
        #pragma omp critical
        #endif
        std::cout << ">> Frame " << frame+1 << " / " << nT
                  << " : Rel. Err. = " << (GainMatrix*v-m).norm()/m.norm()
                  << " : Nb. iter. MinRes = " << niter
                  << std::endl;
    }
}

void MN_inverse (matrice& EstimatedData, const matrice& Data, const matrice& GainMatrix, double SmoothWeight) {
    matrice eye(GainMatrix.nlin(),GainMatrix.nlin());
    eye.set(0);
    for(size_t i = 0; i < GainMatrix.nlin(); ++i) {
        eye(i,i) = SmoothWeight;
    }
    EstimatedData = GainMatrix.transpose() * (GainMatrix * GainMatrix.transpose() + eye).inverse() * Data;
}

// ================= Iterative Mininum norm inversion =======================//

class IMN_inverse_matrice : public virtual matrice
{
public:
    IMN_inverse_matrice (const matrice& Data, const matrice& GainMatrix, double SmoothWeight);
    virtual ~IMN_inverse_matrice () {};
};

IMN_inverse_matrice::IMN_inverse_matrice (const matrice& Data, const matrice& GainMatrix, double SmoothWeight) {
    std::cout << "Running Iterative MN inversion" << std::endl;
    MN_Hessian hess(GainMatrix,SmoothWeight);
    LIN_inverse(*this,hess,GainMatrix,Data);
}

// ================= Mininum norm inversion =======================//

class MN_inverse_matrice : public virtual matrice
{
public:
    MN_inverse_matrice (const matrice& Data, const matrice& GainMatrix, double SmoothWeight);
    virtual ~MN_inverse_matrice () {};
};

MN_inverse_matrice::MN_inverse_matrice (const matrice& Data, const matrice& GainMatrix, double SmoothWeight) {
    std::cout << "Running MN inversion" << std::endl;
    MN_inverse(*this,Data,GainMatrix,SmoothWeight);
}

// ================= Weighted Mininum norm inversion =======================//

class WMN_inverse_matrice : public virtual matrice
{
public:
    WMN_inverse_matrice (const matrice& Data, const matrice& GainMatrix, double SmoothWeight);
    virtual ~WMN_inverse_matrice () {};
};

WMN_inverse_matrice::WMN_inverse_matrice (const matrice& Data, const matrice& GainMatrix, double SmoothWeight) {
    std::cout << "Running WMN inversion" << std::endl;
    WMN_Hessian hess(GainMatrix,SmoothWeight);
    LIN_inverse(*this,hess,GainMatrix,Data);
}

// ================= Gradient based Mininum norm inversion ================ //

class HEAT_inverse_matrice : public virtual matrice
{
public:
    HEAT_inverse_matrice (const matrice& Data, const matrice& GainMatrix, const sparse_matrice& SmoothMatrix, double SmoothWeight);
    virtual ~HEAT_inverse_matrice () {};
};

HEAT_inverse_matrice::HEAT_inverse_matrice (const matrice& Data, const matrice& GainMatrix, const sparse_matrice& SmoothMatrix, double SmoothWeight) {
    std::cout << "Running HEAT inversion" << std::endl;
    fast_sparse_matrice fastSmoothMatrix(SmoothMatrix);
    fast_sparse_matrice fastSmoothMatrix_t(SmoothMatrix.transpose());
    HEAT_Hessian hess(GainMatrix,fastSmoothMatrix,fastSmoothMatrix_t,SmoothWeight);
    LIN_inverse(*this,hess,GainMatrix,Data);
}

// ================= Total variation based inversion =================== //

void TV_inverse(matrice& EstimatedData, const matrice& Data, const matrice& GainMatrix, const sparse_matrice& SmoothMatrix, const vecteur& AiVector, double SmoothWeight, size_t MaxNbIter, double StoppingTol)
{
    fast_sparse_matrice fastSmoothMatrix(SmoothMatrix);
    fast_sparse_matrice fastSmoothMatrix_t(SmoothMatrix.transpose());

    size_t nT = Data.ncol();
    EstimatedData = matrice(GainMatrix.ncol(),nT);

    // #ifdef USE_OMP
    // #pragma omp parallel for
    // #endif
    for(size_t frame=0;frame<nT;frame++)
    {
        cout << ">> Frame " << frame+1 << " / " << nT << endl;
        vecteur m = Data.getcol(frame);
        vecteur v(EstimatedData.nlin());

        // ====================  initialization of source vector ===================== //
        if(frame==0) v.set(0.0);
        else v = EstimatedData.getcol(frame-1);

        double tv_v = compute_tv(v,fastSmoothMatrix,fastSmoothMatrix_t,AiVector);

        bool errorTest = true;

        // ========  Backtracking line search parameters for gradient step  ========= //
        double alpha = 0.001;
        double beta = 0.5;
        int max_iter_line_search = 10;

        // ================== Inverse problem via gradient descent ================== //
        int t;
        for(t=0;t<MaxNbIter && errorTest;t++)
        {
            vecteur gradtv = gentv(v,fastSmoothMatrix,fastSmoothMatrix_t,AiVector);
            vecteur err_vec = GainMatrix*v-m;
            vecteur graddata = GainMatrix.tmult(err_vec);
            vecteur grad = (-SmoothWeight)*gradtv - graddata;
            double f_v = pow(err_vec.norm(),2) + SmoothWeight*tv_v;

            // ======= Backtracking line search for gradient step ======= //
            double search_slope = alpha*grad.norm();
            double f_v_dv;
            double tv_v_dv;
            vecteur v_dv;

            double grad_step = 1.0;

#define USE_LINE_SEARCH 1

#if USE_LINE_SEARCH
            int iter_line_search = 0;
            bool stop_line_search = false;
            while ( stop_line_search != true && (++iter_line_search < max_iter_line_search) ) {
                v_dv = v+grad_step*grad;
                double f_v_dv_data = pow((m-GainMatrix*(v_dv)).norm(),2);
                tv_v_dv = compute_tv(v_dv,fastSmoothMatrix,fastSmoothMatrix_t,AiVector);
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
                       f_v,(err_vec).norm()/m.norm(),tv_v,tol,grad_step,t);
        }
        printf("Total number of iterations : %d\n",t);

        //===========================================================================//
        EstimatedData.setcol(frame,v);
    }
}

class TV_inverse_matrice : public virtual matrice
{
public:
    TV_inverse_matrice (const matrice& Data, const matrice& GainMatrix, const sparse_matrice& SmoothMatrix, const vecteur& AiVector, double SmoothWeight, size_t MaxNbIter, double StoppingTol);
    virtual ~TV_inverse_matrice () {};
};

TV_inverse_matrice::TV_inverse_matrice (const matrice& Data, const matrice& GainMatrix, const sparse_matrice& SmoothMatrix, const vecteur& AiVector, double SmoothWeight, size_t MaxNbIter, double StoppingTol) {
    std::cout << "Running TV inversion" << std::endl;
    TV_inverse(*this,Data,GainMatrix,SmoothMatrix,AiVector,SmoothWeight,MaxNbIter,StoppingTol);
}
