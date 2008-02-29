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

#ifndef H_SYMMATRICE
#define H_SYMMATRICE

#include <sstream>
#include <cmath>
#include <cstdlib>

#include "MatLibConfig.h"
#include "matrice_dcl.h"
#include "symmatrice_dcl.h"
#include "om_utils.h"

inline symmatrice::symmatrice():n(0),t(0),count(0) {}
inline symmatrice::symmatrice(size_t N) { alloc(N); }
inline symmatrice::symmatrice(double* T, int* COUNT, size_t N):n(N),t(T),count(COUNT) {(*count)++;}
inline size_t symmatrice::nlin() const { return n;}
inline size_t symmatrice::ncol() const { return n;}
inline bool symmatrice::empty() const {return t==0;}
inline double* symmatrice::DangerousGetData () {return t;}
inline int* symmatrice::DangerousGetCount () {return count;}
inline double symmatrice::operator()(size_t i,size_t j) const {
    assert(i<n && j<n);
    if(i<=j)
        return t[i+j*(j+1)/2];
    else
        return t[j+i*(i+1)/2];
}

inline double& symmatrice::operator()(size_t i,size_t j) {
    assert(i<n && j<n);
    if(i<=j)
        return t[i+j*(j+1)/2];
    else
        return t[j+i*(i+1)/2];
}

inline void symmatrice::operator /=(double x) {(*this)*=(1/x);}

inline void symmatrice::write(std::ostream& f) const 
{
    f.write((const char*)&n,(std::streamsize)sizeof(int));
    f.write((const char*)t,(std::streamsize)(n*(n+1))/2*sizeof(double));
}

inline void symmatrice::read(std::istream& f) 
{
    destroy();
    f.read((char*)&n,(std::streamsize)sizeof(int));
    alloc(n);
    f.read((char*)t,(std::streamsize)(n*(n+1))/2*sizeof(double));
}

inline void symmatrice::alloc(size_t N)
{
    n=N;
    t=new double[(N*(N+1))/2];
    count=new int[1];
    (*count)=1;
}

inline void symmatrice::destroy()
{
    if (t!=0) {
        (*count)--;
        if ((*count)==0) {
            delete[] t;
            delete[] count;
        }
    }
}

inline void symmatrice::copy(const symmatrice& A)
{
    t=A.t;
    n=A.n;
    if (t) {
        count=A.count;
        (*count)++;
    }
}

inline symmatrice symmatrice::duplicate() const
{
    symmatrice A(n);
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)(n*(n+1))/2,t,1,A.t,1);
#else
    for (size_t i=0;i<(n*(n+1))/2;i++)
        A.t[i]=t[i];
#endif
    return A;
}

inline const symmatrice& symmatrice::operator=(const symmatrice& A) {
    destroy();
    copy(A);
    return *this;
}

inline const symmatrice& symmatrice::operator=(const double d) {
    for(size_t i=0;i<n*n;i++) t[i]=d;
    return *this;
}

inline symmatrice::symmatrice(const symmatrice& A) {
    copy(A);
}

inline symmatrice::symmatrice(const vecteur& v) {
    size_t N = v.size();
    n = (size_t)((sqrt((double)(1+8*N))-1)/2+0.1);
    assert(n*(n+1)/2==N);
    t=v.DangerousGetData();
    if (t) {
        count=v.DangerousGetCount();
        (*count)++;
    }    
}

inline symmatrice::symmatrice(const matrice& A) {
    assert(A.m==A.n);
    alloc(A.n);
    for (size_t j=0; j<n; j++)
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)(j+1),A.t+j*A.m,1,t+(j*(j+1))/2,1);
#else
    for (size_t i=0; i<=j; i++)
        (*this)(i,j)=A(i,j);
#endif
}

inline symmatrice::symmatrice(const char *filename, const char c) {
    t=0;
    if(c=='t') this->loadTxt(filename);
    if(c=='b') this->loadBin(filename);
}

inline void symmatrice::set(double x) {
    for (size_t i=0;i<(n*(n+1))/2;i++)
        t[i]=x;
}

inline vecteur symmatrice::operator *(const vecteur &v) const {
    assert(n==v.n);
    vecteur y(n);
#ifdef HAVE_BLAS
    DSPMV(CblasUpper,(int)n,1.,t,v.t,1,0.,y.t,1);
#else
    for (size_t i=0;i<n;i++) {
        y(i)=0;
        for (size_t j=0;j<n;j++)
            y(i)+=(*this)(i,j)*v(j);
    }
#endif
    return y;
}

inline symmatrice symmatrice::inverse() const {
#ifdef HAVE_LAPACK
    symmatrice invA=duplicate();
    // LU
    int *pivots=new int[n];
    int info;
    DSPTRF('U',invA.n,invA.t,pivots,info);
    // Inverse
    int size=(int)invA.n*64;
    double *work=new double[size];
    DSPTRI('U',invA.n,invA.t,pivots,work,info);

    delete[] pivots;
    delete[] work;
    return invA;
#else
    std::cerr << "!!!!! Inverse not implemented !!!!!" << std::endl;
    exit(1);
#endif
}

inline symmatrice symmatrice::posdefinverse() const {
    // supposes (*this) is definite positive
    symmatrice invA=duplicate();
#ifdef HAVE_LAPACK
    // U'U factorization then inverse
    int info;
    DPPTRF('U',n,invA.t,info);
    DPPTRI('U',n,invA.t,info);
#else
    std::cerr << "Positive definite inverse not defined" << std::endl;
#endif
    return invA;
}

inline double symmatrice::det() {
    symmatrice invA=duplicate();
    double d = 1.0;
#ifdef HAVE_LAPACK
    // Bunch Kaufmqn
    int *pivots=new int[n];
    int info;
    // TUDUtTt
    DSPTRF('U',invA.n,invA.t,pivots,info);
    if(info<0)
        std::cout << "Big problem in det (DSPTRF)" << std::endl;
    for(int i = 0; i<(int) n; i++){
        if(pivots[i] >= 0)
            d *= invA(i,i);
        else // pivots[i] < 0
            if(i < (int) n-1 && pivots[i] == pivots[i+1]){
                d *= (invA(i,i)*invA(i+1,i+1)-invA(i,i+1)*invA(i+1,i));
                i++;
            }
            else
                std::cout << "Big problem in det" << std::endl;
    }
    delete[] pivots;
#else
    std::cerr << "Determinant not defined without LAPACK" << std::endl;
    exit(1);
#endif
    return(d);
}

inline void symmatrice::eigen(matrice & Z, vecteur & D ){
    // performs the complete eigen-decomposition.
    //  (*this) = Z.D.Z'
    // -> eigenvector are columns of the matrix Z.
    // (*this).Z[:,i] = D[i].Z[:,i]
#ifdef HAVE_LAPACK
    symmatrice symtemp = duplicate();
    D = vecteur(n);
    Z = matrice(n,n);

    int info;
    double lworkd;
    int lwork;
    int liwork;

    DSPEVD('V','U',n,symtemp.t,D.t,Z.t,n,&lworkd,-1,&liwork,-1,info);
    lwork = (int) lworkd;
    double * work = new double[lwork];
    int * iwork = new int[liwork];
    DSPEVD('V','U',n,symtemp.t,D.t,Z.t,n,work,lwork,iwork,liwork,info);    

    delete[] work;
    delete[] iwork;
#endif
}

inline matrice symmatrice::operator *(const matrice &B) const {
    assert(n==B.m);
    matrice C(n,B.n);

#ifdef HAVE_BLAS
    matrice D(*this);
    DSYMM(CblasLeft,  CblasUpper
        , (int)n, (int)B.n, 
        1. , D.t, (int)D.n, 
        B.t, (int)B.m, 
        0, C.t,(int)C.m);
#else
    for (size_t j=0;j<B.n;j++) 
        for (size_t i=0;i<n;i++)
        {
            C(i,j)=0;
            for (size_t k=0;k<n;k++)
                C(i,j)+=(*this)(i,k)*B(k,j);
        }
#endif
        return C;
}

inline symmatrice symmatrice::operator +(const symmatrice &B) const {
    assert(n==B.n);
    symmatrice C=duplicate();
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(n*(n+1)/2), 1.0, B.t, 1, C.t , 1);
#else
    for (size_t i=0;i<n*(n+1)/2;i++)
        C.t[i]+=B.t[i];
#endif
    return C;
}

inline symmatrice symmatrice::operator -(const symmatrice &B) const
{
    assert(n==B.n);
    symmatrice C=duplicate();
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(n*(n+1)/2), -1.0, B.t, 1, C.t , 1);
#else
    for (size_t i=0;i<n*(n+1)/2;i++)
        C.t[i]-=B.t[i];
#endif
    return C;
}

inline symmatrice symmatrice::operator *(double x) const {
    symmatrice C(n);
    for (size_t k=0; k<n*(n+1)/2; k++) C.t[k] = t[k]*x;
    return C;
}

inline void symmatrice::operator +=(const symmatrice &B) {
    assert(n==B.n);
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(n*(n+1)/2), 1.0, B.t, 1, t , 1);
#else
    for (size_t i=0;i<n*(n+1)/2;i++)
        t[i]+=B.t[i];
#endif
}

inline void symmatrice::operator *=(double x) {
    for (size_t k=0; k<n*(n+1)/2; k++) t[k] *= x;
}

inline void symmatrice::operator -=(const symmatrice &B) {
    assert(n==B.n);
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(n*(n+1)/2), -1.0, B.t, 1, t , 1);
#else
    for (size_t i=0;i<n*(n+1)/2;i++)
        t[i]+=B.t[i];
#endif
}

inline void symmatrice::saveBin( const char *filename ) const {
    FILE *outfile=fopen(filename,"wb");
    if(outfile==NULL) {std::cout<<"Error opening symmetric matrix binary file: " << filename << std::endl; exit(1);}
    unsigned int ui;
    ui=(unsigned int)n;
    fwrite(&ui,sizeof(unsigned int),1,outfile);
    fwrite(t,sizeof(double),(n*(n+1))/2,outfile);
    fclose(outfile);
}

inline void symmatrice::loadBin( const char *filename ) {
    FILE *infile=fopen(filename,"rb");

    if(infile == NULL) {
        std::cerr<<"Error Opening Matrix File "<<filename<<std::endl;
        exit(1);
    }
    
    unsigned int ui;
    fread(&ui,sizeof(unsigned int),1,infile);
    n=ui;
    if(t!=0) destroy();
    alloc(n);
    fread(t,sizeof(double),(n*(n+1))/2,infile);
    fclose(infile);
}

inline void symmatrice::saveSubBin( const char *filename , size_t i_start, size_t i_end, size_t j_start, size_t j_end) const {
    symmatrice tmpMat(i_end-i_start+1);
    const symmatrice &self=*this;

    for(size_t i=i_start;i<=i_end;i++)
        for(size_t j=j_start;j<=j_end;j++)
        {
            tmpMat(i-i_start,j-j_start)=self(i,j);
        }

        tmpMat.saveBin(filename);
}

inline void symmatrice::saveTxt( const char *filename ) const {
    std::ofstream outfile(filename,std::ios::out);
    if(!outfile.is_open()) { std::cerr<<"Error opening symmetric matrix text file "<<filename<<std::endl;exit(1); }
    for(size_t i=0;i<n;i++)
        for(size_t j=0;j<n;j++)
        {
            outfile<<this->operator ()(i,j);
            if(j!=n-1) outfile<<"\t"; else outfile<<"\n";
        }
}

inline void symmatrice::saveSubTxt( const char *filename , size_t i_start, size_t i_end, size_t j_start, size_t j_end) const {
    std::ofstream outfile(filename,std::ios::out);
    if(!outfile.is_open()) { std::cerr<<"Error opening symmetric matrix text file "<<filename<<std::endl;exit(1); }
    for(size_t i=i_start;i<=i_end;i++)
        for(size_t  j=j_start;j<=j_end;j++)
        {
            outfile<<this->operator ()(i,j);
            if(j!=n-1) outfile<<"\t"; else outfile<<"\n";
        }
}

inline void symmatrice::load( const char *filename ) {
    char extension[128];
    getNameExtension(filename,extension);
    if(!strcmp(extension,"bin") || !strcmp(extension,"BIN")) loadBin(filename);
    else if(!strcmp(extension,"txt") || !strcmp(extension,"TXT")) loadTxt(filename);
    else {
        std::cout << "Warning : Unknown file extension " << std::endl;
        std::cout << "Assuming ASCII format " << std::endl;
        loadTxt(filename);
    }
}

inline void symmatrice::save( const char *filename ) const {
    char extension[128];
    getNameExtension(filename,extension);
    if(!strcmp(extension,"bin") || !strcmp(extension,"BIN")) saveBin(filename);
    else if(!strcmp(extension,"txt") || !strcmp(extension,"TXT")) saveTxt(filename);
    else {
        std::cout << "Warning : Unknown file extension " << std::endl;
        std::cout << "Assuming ASCII format " << std::endl;
        saveTxt(filename);
    }
}

inline void symmatrice::loadTxt( const char *filename ) {
    std::ifstream file(filename);
    if(!file.is_open()) { std::cerr<<"Error reading symmetric matrix text file : "<<filename<<std::endl;exit(1); }
    std::stringstream sst;

    // determine the size of the LineBuffer
    int buf_size=512*1024;
    char *buf=0;
    do
    {
        buf_size*=2;
        if(buf!=0) delete[] buf;
        buf=new char[buf_size];
        file.clear();
        file.seekg(0,std::ios::beg);
        file.getline(buf,buf_size);

    } while(file.fail());

    // Determine the number of columns on the first line
    sst<<buf;
    n=0;
    while(!sst.fail())
    {
        double d; sst>>d;
        n++;
    }
    n--;

    // Determine the number of lines
    file.clear();
    file.seekg(0,std::ios::beg);
    size_t m=0;
    while(!file.fail())
    {
        file.getline(buf,buf_size);
        m++;
    }
    m--;

    assert(m==n);
    if(t!=0) destroy();
    alloc(n);

    // Filling the array
    file.clear();
    file.seekg(0,std::ios::beg);
    for(size_t i=0;i<m;i++)
    {
        std::stringstream s;
        file.getline(buf,buf_size);
        s<<buf;
        for(size_t j=0;j<n;j++)
            s>>(*this)(i,j);
        if(s.fail())
        {
            std::cerr<<"Error Parsing Matrix File "<<filename<<" at line number "<<(int)(i+1)<<std::endl;
            delete[] buf;
            return;
        }
    }
    delete[] buf;
}

inline matrice symmatrice::operator()(size_t i_start, size_t i_end, size_t j_start, size_t j_end) const {
    matrice retMat(i_end-i_start+1,j_end-j_start+1);
    const symmatrice &self=*this;
    for(size_t i=0;i<=i_end-i_start;i++)
        for(size_t j=0;j<=j_end-j_start;j++)
            retMat(i,j)=self(i_start+i,j_start+j);

    return retMat;
}

inline matrice symmatrice::getsubmat(size_t istart, size_t isize, size_t jstart, size_t jsize) const {
    assert ( istart+isize<=n && jstart+jsize<=n );
    return (*this)(istart,istart+isize-1,jstart,jstart+jsize-1);
}

inline symmatrice symmatrice::getsubmat(size_t istart, size_t iend) const {
    assert( iend > istart);
    size_t isize = iend - istart + 1;
    assert ( istart+isize<=n );

    symmatrice mat(isize);
    const symmatrice &self=*this;
    for(size_t i=istart;i<=iend;i++)
        for(size_t j=i;j<=iend;j++)
            mat(i,j)=self(i,j);

    return mat;
}

//returns the solution of (this)*X = B
inline vecteur symmatrice::solveLin(const vecteur &B) const {
    symmatrice invA=duplicate();
    vecteur X = B.duplicate();
    
#ifdef HAVE_LAPACK
    // Bunch Kaufman Factorization
    int *pivots=new int[n];
    int info;
    DSPTRF('U',invA.n,invA.t,pivots,info);
    // Inverse
    int size=(int)invA.n*64;
    double *work=new double[size];
    //dsptrs(uplo,  n    , nrhs , ap   ,ipiv,    b, ldb, info );
    DSPTRS('U',invA.n,1,invA.t,pivots,X.t,invA.n,info);

    delete[] pivots;
    delete[] work;
#else
    std::cout << "solveLin not defined" << std::endl;
#endif
    return X;
}

// stores in B the solution of (this)*X = B, where B is a set of nbvect vector
inline void symmatrice::solveLin(vecteur * B, int nbvect) {
    symmatrice invA=duplicate();
    
#ifdef HAVE_LAPACK
    // Bunch Kaufman Factorization
    int *pivots=new int[n];
    int info;
    //char *uplo="U";
    DSPTRF('U',invA.n,invA.t,pivots,info);
    // Inverse
    int size=(int) invA.n*64;
    double *work=new double[size];
    for(int i = 0; i < nbvect; i++)
        DSPTRS('U',invA.n,1,invA.t,pivots,B[i].t,invA.n,info);

    delete[] pivots;
    delete[] work;
#else
    std::cout << "solveLin not defined" << std::endl;
#endif
}

#endif
