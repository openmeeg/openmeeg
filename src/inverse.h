#define EPSILON 1e-6
#define MINRES_TOL 1e-5

//static double K;
static double k2;
struct tv_evaluator;

inline double tik (const double &x)
{
    return 0.5*x*x;
}
inline double tikp (const double &x)
{
    return x;
}

inline double tikpp (const double &x)
{
    return 1.0; // FIXME:
}

inline double ftv (const double &x)
{
    return x;
}
inline double ftvp (const double &x)
{
    return 1.0; // FIXME:
}

inline double ftvpp (const double &x)
{
    return 0.0; // FIXME:
}

inline double pm (const double &x)
{
    return -0.5*k2*(exp(-(x*x)/k2)-1);
}

inline double pmp (const double &x)
{
    return x*exp(-(x*x)/k2);
}

inline double aub (const double &x)
{
    return (sqrt(1+x*x*k2)-1);
}

inline double aubp (const double &x)
{
    return x/(k2*sqrt((k2+x*x)*k2));
}

inline double aubpp (const double &x)
{
    return 1.0/((k2+x*x)*sqrt((k2+x*x)*k2));
}

inline vecteur gentv( vecteur x,
                      const fast_sparse_matrice &mat,
                      const fast_sparse_matrice &mat_t,
                      const vecteur &Ai, fast_sparse_matrice &hess,
                      double *tv=NULL,
                      double (*f) (const double &)=0,
                      double (*fp) (const double&)=0,
                      double (*fpp) (const double&)=0 )
{ //second order

    // hess is suppposed to be the same size than mat
    vecteur v= mat* x;
    //vecteur zz(v.size());
    vecteur mynorms( v.size()/3);
    vecteur mynorms_inv( v.size()/3);
    for(size_t i=0;i<mynorms.size();i++)
    {
        double *pt=&v(3*i);
        //double *ptz=&zz(3*i);
        mynorms(i)=sqrt(pt[0]*pt[0]+pt[1]*pt[1]+pt[2]*pt[2]);
        mynorms_inv(i)=mynorms(i)!=0?1.0/(mynorms(i)+EPSILON):0;
        //ptz[0]=pt[0]*mynorms_inv(i); ptz[1]=pt[1]*mynorms_inv(i); ptz[2]=pt[2]*mynorms_inv(i);
        double normaliz=mynorms_inv(i);
        pt[0]*=normaliz; pt[1]*=normaliz; pt[2]*=normaliz;
    }

    if(tv!=NULL && f!=0) {*tv=0; for(size_t i=0;i<mynorms.size();i++) *tv+=f(mynorms(i))*Ai(i);}
    if(tv!=NULL && f==0) {*tv=0; for(size_t i=0;i<mynorms.size();i++) *tv+=mynorms(i)*Ai(i);}

    // la hessienne a la même structure creuse que mat c'est à dire 3*v.size() éléments non nuls
    const double* A=mat.getbuf();
    double* z=&v(0);
    double* Hess=hess.getbuf();
    for(size_t k=0;k<v.size()/3;k++)
    {
        double d=mynorms_inv(k);
        //double d2=1.0/(Ai(k)*fp(mynorms(i)));
        Vect3 a1(A[0]*d,A[3]*d,A[6]*d);
        Vect3 a2(A[1]*d,A[4]*d,A[7]*d);
        Vect3 a3(A[2]*d,A[5]*d,A[8]*d);
        //Vect3 zk(z[0]*d2,z[1]*d2,z[2]*d2);
        Vect3 zk(z[0],z[1],z[2]);
        //double nzk=zk.norme();
        double nzk=mynorms(k);
        //double _f= f!=0?f(nzk):nzk;
        double _fp = fp!=0?fp(nzk):1;
        double _fpp = fpp!=0?fpp(nzk)*nzk:0;
        Vect3 alpha1 = ( _fp*a1+( (_fpp-_fp)*(a1*zk) )*zk )*Ai(k);
        Vect3 alpha2 = ( _fp*a2+( (_fpp-_fp)*(a2*zk) )*zk )*Ai(k);
        Vect3 alpha3 = ( _fp*a3+( (_fpp-_fp)*(a3*zk) )*zk )*Ai(k);
        Hess[0]=alpha1[0]; Hess[3]=alpha1[1]; Hess[6]=alpha1[2];
        Hess[1]=alpha2[0]; Hess[4]=alpha2[1]; Hess[7]=alpha2[2];
        Hess[2]=alpha3[0]; Hess[5]=alpha3[1]; Hess[8]=alpha3[2];
        A += 9;
        Hess += 9;
        z += 3;
    }

    for(size_t i=0;i<mynorms.size();i++)
    {
        double *pt=&v(3*i);
        double normaliz=fp!=0?fp(mynorms(i)):1*Ai(i);
        pt[0]*=normaliz; pt[1]*=normaliz; pt[2]*=normaliz;
    }

    return mat_t*v;

}

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


class TvInverseHessian : public LinOp
{
    const matrice &Transfer;
    const fast_sparse_matrice &mat_t;
    const fast_sparse_matrice &HessianAux;
    const double alpha;
    const vecteur *Ai;
public:
    TvInverseHessian(const matrice &TransferMat,
                     const fast_sparse_matrice &mat_t_Hessian,
                     const fast_sparse_matrice &HessianAuxMat,
                     const double &Alpha,
                     const vecteur *AiVec=0):
    Transfer(TransferMat),mat_t(mat_t_Hessian),HessianAux(HessianAuxMat),alpha(Alpha),Ai(AiVec) {}
    virtual vecteur operator * ( const vecteur &x) const
    {
        if(Ai==0)
            return Transfer.tmult(Transfer*x)+alpha*(mat_t*(HessianAux*x));
        else
        {
            vecteur toto=HessianAux*x;
            for(size_t i=0;i<toto.size()/3;i++) {
                toto(3*i)*=(*Ai)(i);
                toto(3*i+1)*=(*Ai)(i);
                toto(3*i+2)*=(*Ai)(i);
            }
            return /*Transfer.tmult(Transfer*x)+*/alpha*(mat_t*toto)/*+1e-6*x*/;
        }
    }
};

class TikInverseHessian : public LinOp
{
    const matrice &Transfer;
    const double alpha;

public:

    TikInverseHessian(const matrice &TransferMat, const double &Alpha):Transfer(TransferMat),alpha(Alpha) {}

    virtual vecteur operator * (const vecteur &x) const
    {
        return Transfer.tmult(Transfer*x)+alpha*x;
    }

};

class GenTvInverseHessian : public LinOp
{
    const matrice &Transfer;
    const fast_sparse_matrice &mat_t;
    const fast_sparse_matrice &HessianAux;
    const double alpha;
    const vecteur *Ai;
public:
    GenTvInverseHessian(const matrice &TransferMat,
        const fast_sparse_matrice &mat_t_Hessian,
        const fast_sparse_matrice &HessianAuxMat,
        const double &Alpha,
        const vecteur *AiVec=0):
    Transfer(TransferMat),mat_t(mat_t_Hessian),HessianAux(HessianAuxMat),alpha(Alpha),Ai(AiVec) {}
    virtual vecteur operator * ( const vecteur &x) const
    {
        if(Ai==0)
            return /*Transfer.tmult(Transfer*x)+alpha**/(mat_t*(HessianAux*x));
        else
        {
            vecteur toto=HessianAux*x;
            for(size_t i=0;i<toto.size()/3;i++)
            {
                toto(3*i)*=(*Ai)(i);
                toto(3*i+1)*=(*Ai)(i);
                toto(3*i+2)*=(*Ai)(i);
            }

            return Transfer.tmult(Transfer*x)+alpha*(mat_t*toto);
        }
    }
};

class Identity : public LinOp
{
public:
    virtual vecteur operator * ( const vecteur &x) const
    {
        return x.duplicate();
    }
};

size_t MinRes2(const LinOp& A,const vecteur& b,vecteur& x0,double tol)
{
    size_t n_max=10000;
    size_t n=1; size_t N=x0.size();
    vecteur v(N); v.set(0.0);
    vecteur v_hat=b-A*x0;
    double beta=v_hat.norm();
    std::cout<<"Residu initial: "<<beta<<std::endl;
    vecteur v_old(v.size());
    vecteur Av(v.size());
    double c=1; double c_old=1; double s_old=0; double s=0;
    vecteur w(N); w.set(0.0);
    vecteur w_oold(N); vecteur w_old=w.duplicate();
    double eta=beta;
    vecteur xMR=x0;//.duplicate();
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
        r1    =sqrt(r1_hat*r1_hat+beta*beta);
        r2    =s_old*alpha+c_oold*c_old*beta_old;
        r3    =s_oold*beta_old;
        //Givens rotation
        c=r1_hat/r1;
        s=beta/r1;
        //update
        w_oold=w_old; w_old=w;
        w=(v-r3*w_oold-r2*w_old)*(1.0/r1);
        xMR+=c*eta*w; norm_rMR=norm_rMR*fabs(s);
        eta=-s*eta;
        // if(!(n%128)) std::cout<<"\r"<<norm_rMR/norm_r0<<std::endl;
    }// end while
    std::cout<<"\r";
    return n;
}

static double (*ftab[4]) (const double &)={0,tik,pm,aub};
static double (*fptab[4]) (const double &)={0,tikp,pmp,aubp};
static double (*fpptab[4]) (const double &)={0,tikpp,0,aubpp};

