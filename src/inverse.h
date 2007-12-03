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
{// first order
    // hess is suppposed to be the same size than mat
    vecteur v= mat* x;
    vecteur mynorms( v.size()/3);
    vecteur mynorms_inv( v.size()/3);
    for(size_t i=0;i<mynorms.size();i++)
    {
        double *pt=&v(3*i);
        mynorms(i)=sqrt(pt[0]*pt[0]+pt[1]*pt[1]+pt[2]*pt[2]);
        mynorms_inv(i)=mynorms(i)!=0?1.0/(mynorms(i)+EPSILON):0;
        double normaliz=mynorms_inv(i)*Ai(i);
        if(fp!=0) normaliz*=fp(mynorms(i));
        pt[0]*=normaliz; pt[1]*=normaliz; pt[2]*=normaliz;
    }

    if(tv!=NULL && f!=0) {*tv=0; for(size_t i=0;i<mynorms.size();i++) *tv+=f(mynorms(i))*Ai(i);}
    if(tv!=NULL && f==0) {*tv=0; for(size_t i=0;i<mynorms.size();i++) *tv+=mynorms(i)*Ai(i);}

    return mat_t*v;
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

        //==========  initialization of source vector =======================//
        vecteur v(GainMatrix.ncol()); v.set(0.0);

        //==========  Invert =======================//
        size_t niter = MinRes2(hess,GainMatrix.tmult(m),v,MINRES_TOL);

        for(size_t i=0;i<EstimatedData.nlin();i++) EstimatedData(i,frame) = v(i);

        #pragma omp critical
        std::cout << ">> Frame " << frame+1 << " / " << nT 
                  << " : Rel. Err. = " << (GainMatrix*v-m).norm()/m.norm() 
                  << " : Nb. iter. MinRes = " << niter
                  << std::endl;
    }
}

//==========  Mininum norm inversion =======================//

class MN_inverse_matrice : public virtual matrice
{
public:
    MN_inverse_matrice (const matrice& Data, const matrice& GainMatrix, double SmoothWeight);
    virtual ~MN_inverse_matrice () {};
};

MN_inverse_matrice::MN_inverse_matrice (const matrice& Data, const matrice& GainMatrix, double SmoothWeight) {
    std::cout << "Running MN inversion" << std::endl;
    MN_Hessian hess(GainMatrix,SmoothWeight);
    LIN_inverse(*this,hess,GainMatrix,Data);
}

//==========  Weighted Mininum norm inversion =======================//

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

//==========  Gradient based Mininum norm inversion =======================//

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

//==========  Total variation based inversion =======================//

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
        else for(size_t i=0;i<v.size();i++) v(i) = EstimatedData(i,frame-1);

        bool errorTest=true;
        double dtv=0.0;

        // ==========================  the inverse problem ========================== //
        int t;
        for(t=0;t<MaxNbIter && errorTest;t++)
        {
            vecteur gradtv = gentv(v,fastSmoothMatrix,fastSmoothMatrix_t,AiVector,&dtv);
            vecteur current_mes = GainMatrix*v;
            vecteur err_vec = current_mes-m;
            vecteur graddata = GainMatrix.tmult(err_vec);
            vecteur Ggraddata = GainMatrix*graddata;

            double denom_data = Ggraddata*Ggraddata;
            double opt_step_data = -(Ggraddata*err_vec)/denom_data;
            vecteur grad = (-SmoothWeight)*gradtv - graddata;
            v = v+grad;
            double tol = sqrt((grad*grad)/(v*v));
            errorTest = tol>StoppingTol;
            if ((t%100)==0)
                printf("TV= %f   Relative Error= %f   opt_step_data= %f   Tol= %f   Iter %d\n",dtv,(err_vec).norm()/m.norm(),opt_step_data,tol,t);
        }
        //===========================================================================//
        for(size_t i=0;i<EstimatedData.nlin();i++) EstimatedData(i,frame)=v(i);

        std::cout << ">> Frame " << frame+1 << " / " << nT  << " : Relative Error = " << (GainMatrix*v-m).norm()/m.norm() << std::endl;
        cout << "\tNumber of iterations = " << t << endl;
        cout << "\tTotal Variation = " << dtv << endl;
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
