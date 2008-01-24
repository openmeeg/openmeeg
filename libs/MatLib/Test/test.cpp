#include "MatLibConfig.h"
#include "vecteur.h"
#include "matrice.h"
#include "symmatrice.h"
#include "sparse_matrice.h"
#include "fast_sparse_matrice.h"
#include "chrono.h"
#include <cmath>
#include <iostream>

using namespace std;

void genericTest( const genericMatrix &M)
{
    cout<<" Generic Test "<<endl;
    cout<<"   nlin  = "<<(int)M.nlin()<<endl;
    cout<<"   ncol  = "<<(int)M.ncol()<<endl;
    cout<<"   (0,0) = "<<M(0,0)<<endl;
    vecteur v(M.ncol());
    v.set(1);
    v=M*v;
    cout<<"   operator * OK"<<endl;
}

int main ()
{
    // section vecteur
    cout<<endl<<"========== vecteurs =========="<<endl;
    vecteur v(8);
    v.set(0);
    v.saveBin("tmp.bin");
    for(int i=0;i<8;i++) v(i)=i;
    v.saveTxt("tmp.txt");
    vecteur vt,vb;
    vb.loadBin("tmp.bin");
    vt.loadTxt("tmp.txt");
    cout<<"v= "<<endl<<v<<endl;
#if defined(USE_ACML)
    cout<<"exp(v)= "<<endl<<v.exp()<<endl;
    cout<<"vb= "<<endl<<vb<<endl;
    cout<<"vt= "<<endl<<vt<<endl;
#endif

    // section matrice
    cout<<endl<<"========== matrices =========="<<endl;
    matrice M(4,4);
    
    for(int i=0;i<4;i++)
        for(int j=0;j<4;j++)
            M(i,j)=pow(2.0,(double)i)+pow(2.0,(double)j);
    
    M.saveBin("tmp.bin");
    M.saveMat("tmp.mat");

    matrice Q=M(1,2,0,2);
    genericTest(Q);

    Q.saveTxt("tmp.txt");
    M.loadTxt("tmp.txt");
    cout<<"Q= "<<endl<<Q<<endl;
    cout<<"M= "<<endl<<M<<endl;
    M.loadBin("tmp.bin");
    cout<<"M= "<<endl<<M<<endl;
    M.loadMat("tmp.mat");
    cout<<"M= "<<endl<<M<<endl;

    ofstream ofs("tmp.test",ios::binary);
    M.write(ofs);
    Q.write(ofs);
    ofs.close();

    matrice tmp;
    ifstream ifs("tmp.test",ios::binary);
    tmp.read(ifs);
    bool flag=true;
    for(size_t i=0;i<M.nlin()*M.ncol();i++) flag=flag && (M[i]==tmp[i]);
    if(flag) cout<<"load/save M is OK"<<endl;
    else cerr<<"load/save M ERROR"<<endl;
    tmp.read(ifs);
    flag=true;
    for(size_t i=0;i<Q.nlin()*Q.ncol();i++) flag=flag && (Q[i]==tmp[i]);
    if(flag) cout<<"load/save Q is OK"<<endl;
    else cerr<<"load/save Q ERROR"<<endl;
    ifs.close();
    
    matrice P(3,3);
    P(0,0) = 25 ; P(0,1) = 3 ; P(0,2) = 6 ;
    P(1,0) = 12 ; P(1,1) = 5 ; P(1,2) = 32 ;
    P(2,0) = 4 ; P(2,1) = 10 ; P(2,2) = 4 ;
    cout << "Matrice P : " << endl;
    for(unsigned int i=0;i<P.nlin(); i++){
      for(unsigned int j=0;j<P.ncol(); j++)
         cout << P(i,j) << "\t" ;
      cout << endl;
    }
    matrice inverseP = P.inverse();
    cout << "Matrice inverse de P : " << endl;
    for(unsigned int i=0;i<inverseP.nlin(); i++){
      for(unsigned int j=0;j<inverseP.ncol(); j++)
         cout << inverseP(i,j) << "\t" ;
      cout << endl;
    }
    matrice unit = P*inverseP;
    double eps = 0.01;
    for(unsigned int i = 0; i<unit.nlin(); i++)
        for(unsigned int j = 0; j<unit.ncol(); j++){
            if(i == j){
                if(abs(unit(i,j)-1) > eps){
                    cerr << "Error: inverse is WRONG-1" << endl;
                    exit(1);
                }
            }
            else
                if(abs(unit(i,j)) > eps){
                    cerr << "Error: inverse is WRONG-2 " << "unit(" << i << "," << j << ") = " << unit(i,j) << endl;
                    exit(1);
                }
        }

    // section symmatrice
    cout<<endl<<"========== symmetric matrices =========="<<endl;
    symmatrice S(4);
    for(unsigned int i=0; i<4; i++)
        for(unsigned int j=i; j<4; j++)
            S(i,j)=pow(2.0,(double)i)+pow(3.0,(double)j);
    
    S.saveBin("tmp.bin");
    genericTest(S);
    matrice R=S(1,2,0,2);
    S.saveTxt("tmp.txt");
    cout<<"S= "<<endl<<S<<endl;
    cout<<"R= "<<endl<<R<<endl;
    R.loadTxt("tmp.txt");
    cout<<"R= "<<endl<<R<<endl;

    // section sparse_matrice 
    cout<<endl<<"========== sparse matrices =========="<<endl; // FIXME : write test on sparse matrices
    sparse_matrice spM(10,10);
    int _n=0;
    for(unsigned int i=0;i<30;i++)
    {
        _n=(_n*1237+1493)%1723;
        int _p=(_n*1237+1493)%1723;
        spM(_n%10,_p%10)=_n;
    }
    spM.refreshNZ();
    genericTest(spM);
    cout<<spM<<endl;
    spM.saveTxt("tmp.txt");
    spM.loadTxt("tmp.txt");
    cout<<spM<<endl;
    spM.saveBin("tmp.bin");
    spM.loadBin("tmp.bin");
    cout<<spM<<endl;
    vecteur spv(10); spv.set(1);
    cout<<"sparse-vector product: "<<spM*spv<<endl;
    cout<<"fast sparse build: "<<endl;
    fast_sparse_matrice fspM(spM);
    genericTest(fspM);
    cout<<"fast_sparse-vector product: "<<fspM*spv<<endl;

// =========
// = SPEED =
// =========
#if 0
    cout<<endl<<"========== speed test =========="<<endl;
    int sm_size = 1000;
    matrice sm(sm_size,sm_size);
    for(size_t j=0;j<sm.ncol();j++)
        for(size_t i=0;i<sm.nlin();i++)
        sm(i,j)=sqrt((double)i)*j*j+sqrt((double)j)*i*i;
    matrice smt=sm.transpose();
    vecteur sv(sm_size);
    for(size_t i=0;i<sv.size();i++)
        sv(i)=i;

    sparse_matrice ssm(10000,10000);
    _n=0;
    for(unsigned int i=0;i<1000;i++)
        for(unsigned int j=0;j<100;j++)
        {
            int _i=rand()%10000;
            int _j=rand()%10000;
            ssm(_i,_j)=sqrt((double)i+j);
        }

    ssm.refreshNZ();
    fast_sparse_matrice fssm(ssm);
    vecteur ssv(10000);
    for(size_t i=0;i<ssv.size();i++)
        ssv(i)=i;

    vecteur resv;
    matrice resm;
    double resd;
    
    timer_start();
    for(int i=0;i<10;i++)
        resm=sm.inverse();
    cout<<"10 LU inversions (size=1000) in "<<timer_lap()<<" seconds"<<endl;

    timer_start();
    for(int i=0;i<1000000;i++)
        resv=sv+sv;
    cout<<"1 000 000 vector additions (size=1000) in "<<timer_lap()<<" seconds"<<endl;

    timer_start();
    for(int i=0;i<1000000;i++)
        resd=sv*sv;
    cout<<"1 000 000 vector dot products (size=1000) in "<<timer_lap()<<" seconds"<<endl;

#if defined(USE_ACML)
    timer_start();
    for(int i=0;i<10000;i++)
        resv=sv.exp();
    cout<<"10 000 vector exp (size=1000) in "<<timer_lap()<<" seconds"<<endl;
#endif

    timer_start();
    for(int i=0;i<1000;i++)
        resv=sm*sv;
    cout<<"1000 matrix-vector products (size=1000) in "<<timer_lap()<<" seconds"<<endl;
    timer_start();
    for(int i=0;i<1000;i++)
        resv=smt.tmult(sv);
    cout<<"1000 matrix-vector cache-aligned products (size=1000) in "<<timer_lap()<<" seconds"<<endl;

    timer_start();
    for(int i=0;i<1000;i++)
        resv=ssm*ssv;
    cout<<"1000 sparse matrix-vector products (size=10.000,nz=100.000) in "<<timer_lap()<<" seconds"<<endl;
    timer_start();
    for(int i=0;i<1000;i++)
        resv=fssm*ssv;
    cout<<"1000 fast sparse matrix-vector products (size=10.000,nz=100.000) in "<<timer_lap()<<" seconds"<<endl;

    timer_start();
    for(int i=0;i<5;i++)
        resm=sm*sm;
    cout<<"5 matrix-matrix products (size=1000) in "<<timer_lap()<<" seconds"<<endl;
    timer_start();
    for(int i=0;i<5;i++)
        resm=smt.tmult(sm);
    cout<<"5 matrix-matrix cache-aligned products (size=1000) in "<<timer_lap()<<" seconds"<<endl;

#endif

// ================
// = BIG_MAT_TEST =
// ================
#if 0
    // 2GB Test
    matrice bigmat(15800,15800);
    vecteur bigvec(15800);
    cout<<"Filling big matrice and big vector..."<<endl;
    for(size_t j=0;j<bigmat.ncol();j++)
        for(size_t i=0;i<bigmat.nlin();i++)
        bigmat(i,j)=sqrt(i)*j*j+sqrt(j)*i*i;
    for(size_t i=0;i<bigvec.size();i++)
        bigvec(i)=i;
    cout<<"Computing 10 big matrix vector products... ";
    timer_start();
    for(size_t i=0;i<10;i++)
        resv=bigmat*bigvec;
    cout<<timer_lap()<<" seconds"<<endl;
    bigmat=matrice();
    bigvec=vecteur();
#endif

// ====================
// = BIG_BIG_MAT_TEST =
// ====================
#if 0
    // 12GB Test
    matrice bigbigmat(38700,38700);
    vecteur bigbigvec(38700);
    cout<<"Filling big big matrice and big big vector..."<<endl;
    for(size_t j=0;j<bigbigmat.ncol();j++)
        for(size_t i=0;i<bigbigmat.nlin();i++)
            bigbigmat(i,j)=sqrt(i)*j*j+sqrt(j)*i*i;
    for(size_t i=0;i<bigbigvec.size();i++)
        bigbigvec(i)=i;
    cout<<"Computing big big matrix vector product... ";
    timer_start();
    resv=bigbigmat*bigbigvec;
    cout<<timer_lap()<<" seconds"<<endl;
    bigbigmat=matrice();
    bigbigvec=vecteur();
#endif

    return 0;
}
