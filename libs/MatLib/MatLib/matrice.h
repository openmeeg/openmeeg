#ifndef H_MATRICE
#define H_MATRICE

#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <cfloat>

#include "matrice_dcl.h"
#include "symmatrice_dcl.h"
#include "vecteur.h"
#include "om_utils.h"

inline matrice::matrice():m(0),n(0),t(0),count(0) {}
inline matrice::matrice(size_t M,size_t N) { alloc(M,N); }
inline matrice::matrice(double* T, int* COUNT, size_t M, size_t N):m(M),n(N),t(T),count(COUNT) {(*count)++;}
inline size_t matrice::nlin() const { return m;}
inline size_t matrice::ncol() const { return n;}
inline bool matrice::empty() const { return t==0;}
inline double* matrice::DangerousGetData () const {return t;}
inline int* matrice::DangerousGetCount () const {return count;}
inline void matrice::DangerousReshape(size_t M, size_t N)
{
    assert(M*N==m*n);
    m=M;
    n=N;
}

inline matrice matrice::duplicate() const
{
    matrice A;
    if (t) {
        A.alloc(m,n);
        copyout(A.t);
    }
    return A;
}

inline void matrice::copyin(const matrice& A)
{
    if (t) {
        assert(m==A.m && n==A.n);
        copyin(A.t);
    }
}
inline matrice matrice::colsref(size_t jstart, size_t jsize) const
{
    assert (jstart+jsize<=n);
    return matrice(t+m*jstart,count,m,jsize);
}
inline vecteur matrice::colref(size_t j) const
{
    assert (j<n);
    return vecteur(t+m*j,count,m);
}
inline vecteur matrice::subcolref(size_t j, size_t istart, size_t isize) const
{
    assert (j<n && istart+isize<=m);
    return vecteur(t+istart+m*j,count,isize);
}
inline double matrice::operator[](size_t i) const
{
    assert(i<m*n);
    return t[i];
}
inline double& matrice::operator[](size_t i)
{
    assert(i<m*n);
    return t[i];
}
inline double matrice::operator()(size_t i,size_t j) const
{
    assert(i<m && j<n);
    return t[i+m*j];
}
inline double& matrice::operator()(size_t i,size_t j)
{
    assert(i<m && j<n);
    return t[i+m*j];
}

inline std::ostream& operator<<(std::ostream& f,const matrice &M) {
    for (size_t i=0;i<M.nlin();i++) {
        for (size_t j=0;j<M.ncol();j++) {
            //f.width(MatriceDisplayElementWidth);
            f << M(i,j) << " ";
        }
        f << std::endl;
    }
    return f;
}

inline void matrice::alloc(size_t M,size_t N)
{
    m = M;
    n = N;
    t = new double[M*N];
    count = new int[1];
    (*count) = 1;
}

inline void matrice::destroy()
{
    if (t!=0) {
        (*count)--;
        if ((*count)==0) {
            delete[] t;
            delete[] count;
        }
    }
}

inline void matrice::copy(const matrice& A)
{
    t=A.t;
    m=A.m;
    n=A.n;
    if (t) {
        count = A.count;
        (*count)++;
    }
}

inline void matrice::copyout(double * p) const {
    if (!t) return;
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)(n*m),t,1,p,1);
#else
    for (size_t i=0;i<m*n;i++)
        p[i]=t[i];
#endif
}

inline void matrice::copyin(const double * p) {
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)(n*m),p,1,t,1);
#else
    for (size_t i=0;i<m*n;i++)
        t[i]=p[i];
#endif
}

inline const matrice& matrice::operator=(const matrice& A)
{
    destroy();
    copy(A);
    return *this;
}

inline const matrice& matrice::set(const double d)
{
    for(size_t i=0;i<n*m;i++) t[i]=d;
    return *this;
}

inline matrice::matrice(const matrice& A)
{
    copy(A);
}

inline matrice::matrice(const symmatrice& A) {
    alloc(A.n,A.n);
    for (size_t j=0; j<n; j++)
        for (size_t i=0; i<m; i++)
            (*this)(i,j)=A(i,j);
}

inline matrice::matrice(const vecteur& v, size_t M, size_t N)
{
    assert(M*N==v.size());
    m=M;
    n=N;
    t=v.DangerousGetData();
    if (t) {
        count=v.DangerousGetCount();
        (*count)++;
    }
}

inline vecteur matrice::operator *(const vecteur &v) const
{
    assert(n==v.n);
    vecteur y(m);
#ifdef HAVE_BLAS
    DGEMV(CblasNoTrans,(int)m,(int)n,1.0,t,(int)m,v.t,1,0.,y.t,1);
#else
    for (size_t i=0;i<m;i++) {
        y(i)=0;
        for (size_t j=0;j<n;j++)
            y(i)+=(*this)(i,j)*v(j);
    }
#endif

    return y;
}

inline matrice matrice::getsubmat(size_t istart, size_t isize, size_t jstart, size_t jsize) const {
    assert (istart+isize<=m && jstart+jsize<=n);
    if (istart==0 && isize==m) return colsref(jstart,jsize).duplicate();

    matrice a(isize,jsize);
    for (size_t j=0; j<jsize; j++)
#ifdef HAVE_BLAS
        BLAS(dcopy,DCOPY)((int)(isize),t+istart+(jstart+j)*m,1,a.t+j*isize,1);
#elif USE_ACML
        dcopy((int)(isize),t+istart+(jstart+j)*m,1,a.t+j*isize,1);
#else
        for (size_t i=0; i<isize; i++)
            a(i,j)=(*this)(istart+i,jstart+j);
#endif
    return a;
}

inline vecteur matrice::getcol(size_t j) const {
    assert(j<n);
    vecteur v(m);
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)m,t+m*j,1,v.t,1);
#else
    for (size_t i=0;i<m;i++) v.t[i]=t[i+m*j];
#endif
    return v;
}

inline vecteur matrice::getlin(size_t i) const {
    assert(i<m);
    vecteur v(n);
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)n,t+i,(int)m,v.t,1);
#else
    for (size_t j=0;j<n;j++) v.t[j]=t[i+m*j];
#endif
    return v;
}

inline void matrice::setcol(size_t j, const vecteur& v) {
    assert(v.size()==m && j<n);
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)m,v.t,1,t+m*j,1);
#else
    for (size_t i=0;i<m;i++) t[i+m*j]=v.t[i];
#endif
}

inline void matrice::setlin(size_t i, const vecteur& v) {
    assert(v.size()==n && i<m);
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)n,v.t,1,t+i,(int)m);
#else
    for (size_t j=0;j<n;j++) t[i+m*j]=v.t[j];
#endif
}

inline vecteur matrice::tmult(const vecteur &v) const
{
    assert(m==v.n);
    vecteur y(n);
#ifdef HAVE_BLAS
    DGEMV(CblasTrans,(int)m,(int)n,1.,t,(int)m,v.t,1,0.,y.t,1);
#else
    for (size_t i=0;i<n;i++) {
        y(i)=0;
        for (size_t j=0;j<m;j++)
            y(i)+=(*this)(j,i)*v(j);
    }
#endif

    return y;
}

inline matrice matrice::inverse() const
{
#ifdef HAVE_LAPACK
    assert(m==n);
    matrice invA=duplicate();
    // LU
    int *pivots=new int[n];
    int info;
    DGETRF(invA.m,invA.n,invA.t,invA.m,pivots,info);
    // Inverse
    int size=(int)invA.n*64;
    double *work=new double[size];
    DGETRI(invA.n,invA.t,invA.n,pivots,work,size,info);
    delete[] pivots;
    delete[] work;
    return invA;
#else
    std::cerr << "!!!!! Inverse not implemented !!!!!" << std::endl;
    exit(1);
#endif
}

inline matrice matrice::pinverse(double tolrel) const {
#if defined(HAVE_BLAS) && defined(HAVE_LAPACK)
    if(n > m) return transpose().pinverse().transpose();
    else {
        matrice result(n,m);
        matrice U,S,V;
        svd(U,S,V);
        double maxs=0;
        int mimi=(int)std::min(S.nlin(),S.ncol());
        for(int i=0;i<mimi;i++) maxs=std::max(S(i,i),maxs);
        if (tolrel==0) tolrel=DBL_EPSILON;
        double tol = std::max(m,n) * maxs * tolrel;
        int r=0; for(int i=0;i<mimi;i++) if(S(i,i)>tol) r++;
        if (r == 0)
        {
            result.set(0.);
            return result;
        }
        else
        {
            matrice s(r,r); s.set(0);
            for(int i=0;i<r;i++) s(i,i)=1.0/S(i,i);
            matrice Vbis;
            Vbis.DangerousBuild(V.t,V.nlin(),r);
            matrice Ubis;
            Ubis.DangerousBuild(U.t,U.nlin(),r);
            result=Vbis*s*Ubis.transpose();
            Vbis.DangerousKill();
            Ubis.DangerousKill();
            return result;
        }
    }
#else
    std::cerr << "pinv not implemented without blas/lapack" << std::endl;
    exit(1);
#endif
}

inline matrice matrice::transpose() const {
    matrice result(n,m);
    for(size_t i=0;i<m;i++) for(size_t j=0;j<n;j++) result(j,i)=(*this)(i,j);
    return result;
}

inline void matrice::svd(matrice &U,matrice &S, matrice &V) const {
#ifdef HAVE_LAPACK
    matrice cpy=duplicate();
    int mimi=(int)std::min(m,n);
    U=matrice(m,n); U.set(0);
    V=matrice(n,n); V.set(0);
    S=matrice(n,n); S.set(0);
    double *s=new double[mimi];
    int lwork=4 *mimi*mimi + (int)std::max( m, n) + 9*mimi;
    double *work=new double[lwork];
    int *iwork=new int[8*mimi];
    int info;
    DGESDD('S',m,n,cpy.t,m,s,U.t,U.m,V.t,V.m,work,lwork,iwork,info);
    for(int i=0;i<mimi;i++) S(i,i)=s[i];
    V=V.transpose();
    delete[] s;
    delete[] work;
    delete[] iwork;
#else
    std::cerr<<"svd not implemented without blas/lapack"<<std::endl;
#endif
}

inline matrice matrice::operator *(const matrice &B) const
{
    assert(n==B.m);
    size_t p=n;
    matrice C(m,B.n);
#ifdef HAVE_BLAS
    DGEMM(CblasNoTrans,CblasNoTrans,
        (int)C.m,(int)C.n,(int)p,
        1.,t,(int)m,
        B.t,(int)B.m,
        0.,C.t,(int)C.m);
#else
    for (size_t i=0;i<C.m;i++)
        for (size_t j=0;j<C.n;j++) {
            C(i,j)=0;
            for (size_t k=0;k<p;k++)
                C(i,j)+=(*this)(i,k)*B(k,j);
        }
#endif
        return C;
}

inline matrice matrice::tmult(const matrice &B) const
{
    assert(m==B.m);
    size_t p=m;
    matrice C(n,B.n);
#ifdef HAVE_BLAS
    DGEMM(CblasTrans,CblasNoTrans,
        (int)C.m,(int)C.n,(int)p,
        1.,t,(int)m,
        B.t,(int)B.m,
        0.,C.t,(int)C.m);
#else
    for (size_t i=0;i<C.m;i++)
        for (size_t j=0;j<C.n;j++) {
            C(i,j)=0;
            for (size_t k=0;k<p;k++)
                C(i,j)+=(*this)(k,i)*B(k,j);
        }
#endif
        return C;
}
inline matrice matrice::multt(const matrice &B) const
{
    assert(n==B.n);
    size_t p=n;
    matrice C(m,B.m);
#ifdef HAVE_BLAS
    DGEMM(CblasNoTrans,CblasTrans,
        (int)C.m,(int)C.n,(int)p,
        1.,t,(int)m,
        B.t,(int)B.m,
        0.,C.t,(int)C.m);
#else
    for (size_t i=0;i<C.m;i++)
        for (size_t j=0;j<C.n;j++) {
            C(i,j)=0;
            for (size_t k=0;k<p;k++)
                C(i,j)+=(*this)(i,k)*B(j,k);
        }
#endif
        return C;
}

inline matrice matrice::tmultt(const matrice &B) const
{
    assert(m==B.n);
    size_t p=m;
    matrice C(n,B.m);
#ifdef HAVE_BLAS
    DGEMM(CblasTrans,CblasTrans,
        (int)C.m,(int)C.n,(int)p,
        1.,t,(int)m,
        B.t,(int)B.m,
        0.,C.t,(int)C.m);
#else
    for (size_t i=0;i<C.m;i++)
        for (size_t j=0;j<C.n;j++) {
            C(i,j)=0;
            for (size_t k=0;k<p;k++)
                C(i,j)+=(*this)(k,i)*B(j,k);
        }
#endif
        return C;
}


inline matrice matrice::operator *(const symmatrice &B) const
{
    assert(n==B.n);
    matrice C(m,B.n);

#ifdef HAVE_BLAS
    matrice D(B);
    DSYMM(CblasRight,  CblasUpper
        , (int)m, (int)D.n,
        1. , D.t, (int)D.n,
        t, (int)m,
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

inline matrice matrice::operator *(double x) const {
    matrice C(m,n);
    for (size_t k=0; k<m*n; k++) C.t[k] = t[k]*x;
    return C;
}
inline matrice matrice::operator /(double x) const {
    matrice C(m,n);
    for (size_t k=0; k<m*n; k++) C.t[k] = t[k]/x;
    return C;
}

inline matrice matrice::operator +(const matrice &B) const
{
    assert(n==B.n);
    assert(m==B.m);
    matrice C=duplicate();
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(n*m), 1.0, B.t, 1, C.t , 1);
#else
    for (size_t i=0;i<m*n;i++)
        C.t[i]+=B.t[i];
#endif
    return C;
}

inline matrice matrice::operator -(const matrice &B) const
{
    assert(n==B.n);
    assert(m==B.m);
    matrice C=duplicate();
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(n*m), -1.0, B.t, 1, C.t , 1);
#else
    for (size_t i=0;i<m*n;i++)
        C.t[i]-=B.t[i];
#endif
    return C;
}

inline void matrice::operator +=(const matrice &B)
{
    assert(n==B.n);
    assert(m==B.m);
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(n*m), 1.0, B.t, 1, t , 1);
#else
    for (size_t i=0;i<m*n;i++)
        t[i]+=B.t[i];
#endif
}

inline void matrice::operator -=(const matrice &B)
{
    assert(n==B.n);
    assert(m==B.m);
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)(n*m), -1.0, B.t, 1, t , 1);
#else
    for (size_t i=0;i<m*n;i++)
        t[i]-=B.t[i];
#endif
}

inline void matrice::operator *=(double x) {
    for (size_t k=0; k<m*n; k++) t[k] *= x;
}
inline void matrice::operator /=(double x) {
    for (size_t k=0; k<m*n; k++) t[k] /= x;
}

inline void matrice::saveBin( const char *filename ) const
{
    FILE *outfile=fopen(filename,"wb");
    if(outfile==NULL) {std::cout<<"Error opening matrix binary file: " << filename << std::endl; exit(1);}
    unsigned int ui;
    ui=(unsigned int)m;
    fwrite(&ui,sizeof(unsigned int),1,outfile);
    ui=(unsigned int)n;
    fwrite(&ui,sizeof(unsigned int),1,outfile);
    fwrite(t,sizeof(double),n*m,outfile);
    fclose(outfile);
}

inline void matrice::loadBin( const char *filename )
{
    FILE *infile=fopen(filename,"rb");

    if(infile == NULL) {
        std::cerr<<"Error Opening Matrix File "<<filename<<std::endl;
        exit(1);
    }


    unsigned int ui;
    fread(&ui,sizeof(unsigned int),1,infile);
    m=ui;
    fread(&ui,sizeof(unsigned int),1,infile);
    n=ui;
    if(t!=0) destroy();
    alloc(m,n);
    fread(t,sizeof(double),n*m,infile);
    fclose(infile);
}

inline void matrice::saveSubBin( const char *filename, size_t i_start, size_t i_end, size_t j_start, size_t j_end) const
{
    FILE *outfile=fopen(filename,"wb");
    if(outfile==NULL) {std::cout<<"Error opening matrix binary file: " << filename << std::endl; exit(1);}
    size_t sub_m=i_end-i_start+1;
    size_t sub_n=j_end-j_start+1;

    unsigned int ui;
    ui=(unsigned int)sub_m;
    fwrite(&ui,sizeof(unsigned int),1,outfile);
    ui=(unsigned int)sub_n;
    fwrite(&ui,sizeof(unsigned int),1,outfile);
    for(size_t j=j_start;j<=j_end;j++)
        fwrite(t+i_start+j*m,sizeof(double),sub_m,outfile);

    fclose(outfile);
}

inline void matrice::saveTxt( const char *filename ) const
{
    std::ofstream outfile(filename);
    if(!outfile.is_open()) {std::cerr<<"Error opening matrix text file: " << filename << std::endl; exit(1);}
    for(size_t i=0;i<m;i++)
        for(size_t j=0;j<n;j++)
        {
            outfile<<this->operator ()(i,j);
            if(j!=n-1) outfile<<"\t"; else outfile<<"\n";
        }
}

inline void matrice::saveSubTxt( const char *filename, size_t i_start, size_t i_end, size_t j_start, size_t j_end) const
{
    std::ofstream outfile(filename,std::ios::out);
    if(!outfile.is_open()) {std::cerr<<"Error opening matrix text file: " << filename << std::endl; exit(1);}
    for(size_t i=i_start;i<=i_end;i++)
        for(size_t j=j_start;j<=j_end;j++)
        {
            outfile<<this->operator ()(i,j);
            if(j!=n-1) outfile<<"\t"; else outfile<<"\n";
        }

}

inline void matrice::loadMat(const char *filename) throw(std::string)
{
#ifdef USE_MATIO
    mat_t* mat = Mat_Open(filename,MAT_ACC_RDONLY);
    if (mat) {
        matvar_t* matvar = Mat_VarReadNext(mat);
        while (matvar!=NULL && (matvar->rank!=2 && matvar->data_type!=MAT_T_DOUBLE))
            matvar = Mat_VarReadNext(mat);
        if (matvar==NULL)
            throw std::string("There is no 2D double matrix in file ")+filename;
        m = matvar->dims[0];
        n = matvar->dims[1];
        t = static_cast<double*>(matvar->data);
        matvar->mem_conserve = 1;
        Mat_VarFree(matvar);
        Mat_Close(mat);
    }
#else
    std::cerr << "You have to compile OpenMEEG with MATIO to load matlab files" << std::endl;
#endif
}

inline void matrice::saveMat( const char *filename ) const
{
#ifdef USE_MATIO
    mat_t* mat = Mat_Open(filename,MAT_ACC_RDWR);
    if (mat) {
        matvar_t* matvar;
        int dims[2] = { m, n };
        matvar = Mat_VarCreate("",MAT_C_DOUBLE,MAT_T_DOUBLE,2,dims,t,0);
        Mat_VarWrite(mat,matvar,COMPRESSION_ZLIB);
        Mat_VarFree(matvar);
        Mat_Close(mat);
    }

#else
    std::cerr << "You have to compile OpenMEEG with MATIO to save matlab files" << std::endl;
#endif
}

inline matrice matrice::operator()(size_t i_start, size_t i_end, size_t j_start, size_t j_end) const
{
    matrice retMat(i_end-i_start+1,j_end-j_start+1);

    for(size_t j=0;j<=j_end-j_start;j++)
#ifdef HAVE_BLAS
        BLAS(dcopy,DCOPY)((int)(i_end-i_start+1), t+(j*m)+i_start, 1, retMat.t+(j*retMat.m), 1);
#endif

    return retMat;
}

inline void matrice::load( const char *filename ) {
    char extension[128];
    getNameExtension(filename,extension);
    if(!strcmp(extension,"bin") || !strcmp(extension,"BIN")) loadBin(filename);
    else if(!strcmp(extension,"txt") || !strcmp(extension,"TXT")) loadTxt(filename);
    else if(!strcmp(extension,"mat") || !strcmp(extension,"MAT")) loadMat(filename);
    else {
        std::cout << "Warning : Unknown file extension " << std::endl;
        std::cout << "Assuming ASCII format " << std::endl;
        loadTxt(filename);
    }
}

inline void matrice::loadTxt( const char *filename )
{
    std::ifstream file(filename);
    if(!file.is_open()) {
        std::cerr<<"Error Opening Matrix File "<<filename<<std::endl;
        exit(1);
    }
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
    m=0;
    while(!file.fail())
    {
        file.getline(buf,buf_size);
        m++;
    }
    m--;

    if(t!=0) destroy();
    alloc(m,n);

    // Filling the array
    file.clear();
    file.seekg(0,std::ios::beg);
    for(size_t i=0;i<m;i++)
    {
        std::stringstream s;
        file.getline(buf,buf_size);
        s<<buf;
        for(size_t j=0;j<n;j++)
            s>>t[i+j*m];
        if(s.fail())
        {
            std::cerr<<"Error Parsing Matrix File "<<filename<<" at line number "<<(int)i+1<<std::endl;
            delete[] buf;
            return;
        }
    }

    delete[] buf;
}

inline matrice::matrice(const char* filename, const char c)
{
    t=0;
    if(c=='t') this->loadTxt(filename);
    if(c=='b') this->loadBin(filename);
}

inline void matrice::DangerousBuild( double *pt, size_t i, size_t j)
{
    t=pt;
    m=i;
    n=j;
    count=new int[1];
    *count=1;
}
inline void matrice::DangerousKill ()
{
    delete[] count;
    t=0;
}

inline matrice matrice::lowercholesky() const {
    assert (n==m);
    matrice b=duplicate();
#ifdef HAVE_LAPACK
    int info;
    DPOTF2('L',n,b.t,n,info);
    // il faut mettre encore les zeros a la main
    // dans la partie sup. droite qui n'a pas ÈtÈ overwritten
    for (size_t j=1; j<n; j++) for (size_t i=0; i<j; i++) b(i,j)=0;
#else
    std::cerr << "lowercholesky not defined" << std::endl;
#endif

    return b;
}

inline double matrice::dot(const matrice& b) const {
    assert(nlin()==b.nlin()&&ncol()==b.ncol());
#ifdef HAVE_BLAS
    return BLAS(ddot,DDOT)((int)(nlin()*ncol()),t,1,b.t,1);
#else
    double s=0;
    for (size_t i=0;i<nlin()*ncol();i++)
        s+=t[i]*b.t[i];
    return s;
#endif
}

inline double matrice::frobenius_norm() const {
#ifdef HAVE_LAPACK
    double info;
    matrice b=duplicate();
    return DLANGE('F',m,n,b.t,m,&info);
#else
    double d=0;
    for (size_t i=0; i<m*n; i++) d+=t[i]*t[i];
    return sqrt(d);
#endif
}

inline double matrice::det() const{
    assert (n==ncol());
#ifdef HAVE_LAPACK
    matrice b=duplicate();
    int info;
    // using cholesky factorization
    DPOTF2('L',n,b.t,n,info);
    double d=1;
    for (size_t i=1; i<n; i++) d*=b(i,i);
    return d;
#else
    std::cerr << "lowercholesky not defined" << std::endl;
    return 0;
#endif
}

inline vecteur matrice::mean() const {
    vecteur v(ncol()); v.set(0);
    for(size_t j = 0; j < ncol(); ++j) {
        for(size_t i = 0; i < nlin(); ++i) {
            v(j) += this->operator()(i,j);
        }
    }
    for(size_t j = 0; j < ncol(); ++j) {
        v(j) = v(j) / nlin();
    }
    return v;
}

inline vecteur matrice::tmean() const {
    vecteur v(nlin()); v.set(0);
    for(size_t j = 0; j < ncol(); ++j) {
        for(size_t i = 0; i < nlin(); ++i) {
            v(i) += this->operator()(i,j);
        }
    }
    for(size_t i = 0; i < nlin(); ++i) {
        v(i) = v(i) / ncol();
    }
    return v;
}

inline void matrice::info() const {
    if ((nlin() == 0) && (ncol() == 0)) {
        std::cout << "Matrix Empty" << std::endl;
        return;
    }

    std::cout << "Dimensions : " << nlin() << " x " << ncol() << std::endl;

    double minv = this->operator()(0,0);
    double maxv = this->operator()(0,0);
    size_t mini = 0;
    size_t maxi = 0;
    size_t minj = 0;
    size_t maxj = 0;

    for(size_t i = 0; i < nlin(); ++i)
    {
        for(size_t j = 0; j < ncol(); ++j)
        {
            if (minv > this->operator()(i,j)) {
                minv = this->operator()(i,j);
                mini = i;
                minj = j;
            } else if (maxv < this->operator()(i,j)) {
                maxv = this->operator()(i,j);
                maxi = i;
                maxj = j;
            }
        }
    }
    std::cout << "Min Value : " << minv << " (" << mini << "," << minj << ")" << std::endl;
    std::cout << "Max Value : " << maxv << " (" << maxi << "," << maxj << ")" << std::endl;
    std::cout << "First Values" << std::endl;
    for(size_t i = 0; i < std::min(nlin(),(size_t) 5); ++i)
    {
        for(size_t j = 0; j < std::min(ncol(),(size_t) 5); ++j)
        {
            std::cout << this->operator()(i,j) << " " ;
        }
        std::cout << std::endl ;
    }
}

#endif
