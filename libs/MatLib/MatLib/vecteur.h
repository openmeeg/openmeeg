#ifndef H_VECTEUR
#define H_VECTEUR

#include "MatLibConfig.h"
#include "vecteur_dcl.h"
#include "matrice.h"

inline vecteur::vecteur():n(0),t(0),count(0) {}
inline vecteur::vecteur(size_t N) {alloc(N); }
inline vecteur::vecteur(double* T, int* COUNT, size_t N):n(N),t(T),count(COUNT) {(*count)++;}
inline size_t vecteur::size() const { return n;}
inline bool vecteur::empty() const {return t==0;}
inline double*  vecteur::DangerousGetData () const {return t;}
inline int* vecteur::DangerousGetCount () const {return count;}
inline vecteur vecteur::duplicate() const
{
    vecteur A;
    if (t)
    {
        A.alloc(n);
        copyout(A.t);
    }
    return A;
}
inline void vecteur::copyin(const vecteur& A)
{
    assert(n==A.n);
    copyin(A.t);
}
inline vecteur vecteur::subvectref(size_t istart, size_t s) const {
    assert(istart+s<=n);
    return vecteur(t+istart,count,s);
}
inline double vecteur::operator()(size_t i) const {
    assert(i<n);
    return t[i];
}
inline double& vecteur::operator()(size_t i) {
    assert(i<n);
    return t[i];
}
inline void vecteur::operator/=(double x) {(*this)*=(1/x);}
inline vecteur vecteur::operator/(double x) const {return (*this)*(1/x);}
inline double vecteur::mean() const {return sum()/size();}
inline void vecteur::write(std::ostream& f) const
{
    f.write((const char*)t,(std::streamsize)n*sizeof(double));
}
inline void vecteur::read(std::istream& f)
{
    f.read((char*)t,(std::streamsize)n*sizeof(double));
}
inline vecteur operator * (const double &d, const vecteur &v) {return v*d;}

inline std::ostream& operator<<(std::ostream& f,const vecteur &M) {
    for (size_t i=0;i<M.size();i++)
    {
        f << M(i) << " ";
    }
    return f;
}

inline std::istream& operator>>(std::istream& f,vecteur &M) {
    for (size_t i=0;i<M.size();i++)
    {
        f >> M(i);
    }

    return f;
}

inline void vecteur::alloc(size_t N) {
    n=N;
    t=new double[N];
    count=new int[1];
    (*count)=1;
}

inline void vecteur::destroy() {
    if (t!=0)
    {
        (*count)--;
        if (*count == 0)
        {
            delete[] t;
            delete[] count;
        }
    }
}

inline void vecteur::copy(const vecteur& A) {
    t=A.t;
    n=A.n;
    if (t)
    {
        count=A.count;
        (*count)++;
    }
}

inline void vecteur::copyout(double * p) const {
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)n,t,1,p,1);
#else
    for( size_t i=0; i<n; i++ )
    {
        p[i]=t[i];
    }
#endif
}

inline void vecteur::copyin(const double * p) {
#ifdef HAVE_BLAS
    BLAS(dcopy,DCOPY)((int)n,p,1,t,1);
#else
    for( size_t i=0; i<n; i++ )
    {
        t[i]=p[i];
    }
#endif
}

inline const vecteur& vecteur::operator=(const vecteur& A)
{
    destroy();
    copy(A);
    return *this;
}

inline vecteur::vecteur(const vecteur& A)
{
    copy(A);
}

inline vecteur::vecteur(matrice& A)
{
    n=A.nlin()*A.ncol();
    t=A.DangerousGetData();
    if (t)
    {
        count=A.DangerousGetCount();
        (*count)++;
    }
}

inline vecteur::vecteur(symmatrice& A)
{
    n=A.nlin()*(A.nlin()+1)/2;
    t=A.DangerousGetData();
    if (t)
    {
        count=A.DangerousGetCount();
        (*count)++;
    }
}

inline void vecteur::saveTxt(const char* filename) const
{
    std::ofstream outfile(filename,std::ios::out);
    if(!outfile.is_open()) { std::cerr<<"Error opening vector text file "<<filename<<std::endl;exit(1); }
    for(size_t i=0;i<n;i++) outfile<<t[i]<<std::endl;
}

inline void vecteur::saveBin(const char* filename) const
{
    FILE *f=fopen(filename,"wb");
    if(f==NULL) {std::cerr<<"Error opening vector binary file: " << filename << std::endl; exit(1);}
    unsigned int ui;
    ui=(unsigned int)n;
    fwrite(&ui,sizeof(unsigned int),1,f);
    fwrite(t,sizeof(double),n,f);
    fclose(f);
}

inline void vecteur::loadTxt(const char* filename)
{
    std::ifstream infile(filename,std::ios::in);
    if(!infile.is_open()) {std::cerr<<"Error reading vector text file"<<filename<<std::endl; exit(1);}

    const int buf_size=1024;
    char buf[buf_size];

    //Let's count the number of lines
    n=0;
    while ( (! infile.eof()) ) {infile.getline (buf,buf_size); n++;}
    n--;

    if(t!=0) destroy();
    alloc(n);

    infile.clear();
    infile.seekg(0);
    for(size_t i=0;i<n;i++) infile>>t[i];
}

inline void vecteur::loadBin(const char* filename)
{
    FILE *f=fopen(filename,"rb");
    if(f==NULL) {std::cerr<<"Error reading vector binary file "<<filename<<std::endl; exit(1);}
    unsigned int ui;
    fread(&ui,sizeof(unsigned int),1,f);
    n=ui;
    if(t!=0) destroy();
    alloc(n);
    fread(t,sizeof(double),n,f);
    fclose(f);
}

inline vecteur vecteur::operator+(const vecteur& v) const {
    assert(n==v.n);
    vecteur p=duplicate();
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)n,1,v.t,1,p.t,1);
#else
    for( size_t i=0; i<n; i++ )
    {
        p.t[i]=t[i]+v.t[i];
    }
#endif
    return p;
}

inline vecteur vecteur::operator-(const vecteur& v) const {
    assert(n==v.n);
    vecteur p=duplicate();
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)n,-1,v.t,1,p.t,1);
#else
    for( size_t i=0; i<n; i++ )
    {
        p.t[i]=t[i]-v.t[i];
    }
#endif
    return p;
}

inline void vecteur::operator+=(const vecteur& v) {
    assert(n==v.n);
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)n,1,v.t,1,t,1);
#else
    for( size_t i=0; i<n; i++ )
    {
        t[i]+=v.t[i];
    }

#endif
}

inline void vecteur::operator-=(const vecteur& v) {
    assert(n==v.n);
#ifdef HAVE_BLAS
    BLAS(daxpy,DAXPY)((int)n,-1,v.t,1,t,1);
#else
    for( size_t i=0; i<n; i++ )
    {
        t[i]-=v.t[i];
    }
#endif
}

inline double vecteur::operator*(const vecteur& v) const {
    assert(n==v.n);
#ifdef HAVE_BLAS
    return BLAS(ddot,DDOT)((int)n,t,1,v.t,1);
#else
    double s=0;
    for( size_t i=0; i<n; i++ )
    {
        s+=t[i]*v.t[i];
    }
    return s;
#endif
}

inline vecteur vecteur::kmult(const vecteur& v) const { // Kronecker multiplication
    assert(n == v.n);
// #ifdef HAVE_BLAS
//     // FIXME : add blas version
// #else
    vecteur p(n);
    for( size_t i=0; i<n; i++ )
    {
        p(i) = v(i)*t[i];
    }
// #endif
    return p;
}

inline vecteur vecteur::operator*(double x) const {
#ifdef HAVE_BLAS
    vecteur p=duplicate();
    BLAS(dscal,DSCAL)((int)n,x,p.t,1);
#else
    vecteur p(n);
    for( size_t i=0; i<n; i++ )
    {
        p.t[i]=x*t[i];
    }
#endif
    return p;
}

inline vecteur vecteur::operator+(double x) const
{
    vecteur p=duplicate();
    for( size_t i=0; i<n; i++ )
    {
        p.t[i]+=x;
    }
    return p;
}

inline vecteur vecteur::operator-(double x) const
{
    vecteur p=duplicate();
    for( size_t i=0; i<n; i++ )
        p.t[i]-=x;

    return p;
}

inline void vecteur::operator*=(double x) {
#ifdef HAVE_BLAS
    BLAS(dscal,DSCAL)((int)n,x,t,1);
#else
    for( size_t i=0; i<n; i++ )
    {
        t[i]*=x;
    }
#endif
}

inline double vecteur::norm() const
{
#ifdef HAVE_BLAS
    return BLAS(dnrm2,DNRM2)((int)n,t,1);
#else
    std::cout << "'vecteur::norm' not implemented" << std::endl;
    exit(1);
    return 0;
#endif
}

inline void vecteur::set(double x) {
    assert(n>0);
    for( size_t i=0; i<n; i++ )
    {
        t[i]=x;
    }
}

inline double vecteur::sum() const
{
    double s=0;
    for (size_t i=0; i<n; i++)
    {
        s+=t[i];
    }
    return s;
}

inline void vecteur::DangerousBuild(double *ptr, size_t s)
{
    t=ptr;
    n=s;
    count=new int[1];
    *count=1;
}

inline void vecteur::DangerousKill()
{
    t=0;
    n=0;
    delete[] count;
}


inline vecteur vecteur::conv(const vecteur& v) const {
    if (v.n<n) return v.conv(*this);

    vecteur p(n+v.n-1);
    p.set(0);
    for (size_t i=0; i<v.n; i++) {
#ifdef HAVE_BLAS
        BLAS(daxpy,DAXPY)((int)n, v(i), t, 1, p.t+i, 1);
#else
        for (size_t j=0;j<n;j++)
        {
            p(i+j)+=v(i)*t[j];
        }
#endif
    }
    return p;
}

inline vecteur vecteur::conv_trunc(const vecteur& v) const {
    vecteur p(v.n);
    p.set(0);
    for (size_t i=0; i<v.n; i++)
    {
        size_t m = std::min(n,v.n-i);
#ifdef HAVE_BLAS
        BLAS(daxpy,DAXPY)((int)m, v(i), t, 1, p.t+i, 1);
#else
        for (size_t j=0;j<m;j++)
        {
            p(i+j)+=v(i)*t[j];
        }
#endif
    }
    return p;
}


inline matrice vecteur::outer_product(const vecteur& v) const {
    assert(n==v.size());
    matrice A(n,n);
#ifdef HAVE_BLAS
    DGEMM(CblasNoTrans,CblasNoTrans,
        (int)n,(int)n,1,
        1.,t,(int)n,
        v.t,(int)n,
        0.,A.DangerousGetData(),(int)n);
#else
    for( unsigned int j=0; j<n; j++ ) {
        for ( unsigned int i=0; i<n; i++) {
            A(i,j)=v(i)*(*this)(j);
        }
    }
#endif
    return A;
}

#endif
