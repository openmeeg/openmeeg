#ifndef H_MATRICE_DCL
#define H_MATRICE_DCL

#include "MatLibConfig.h"
#include <iostream>
#include <fstream>
#include <cassert>

#include "generic_matrix.h"

class symmatrice;
class vecteur;

class matrice : public virtual genericMatrix {
protected:
    size_t m,n;
    double *t;
    int* count;

    inline void alloc(size_t M,size_t N);
    inline void copy(const matrice& A);
    inline void destroy();
    inline void copyout(double * p) const;
    inline void copyin(const double * p);
public:
    inline matrice();
    inline matrice(size_t M,size_t N);
    inline matrice(const matrice& A);
    inline explicit matrice(const symmatrice& A);
    inline matrice(const vecteur& v, size_t M, size_t N); // accede violemment a v.t et v.count
    inline matrice(double* T, int* COUNT, size_t M, size_t N);
    inline matrice(const char *filename, char c='t');

    inline ~matrice() { destroy(); }
    inline size_t nlin() const;
    inline size_t ncol() const;
    inline bool empty() const;
    inline void DangerousBuild( double *, size_t i, size_t j);
    inline void DangerousKill ();
    inline double* DangerousGetData () const;
    inline int* DangerousGetCount () const;
    inline void DangerousReshape(size_t M, size_t N);

    inline matrice duplicate() const;
    inline void copyin(const matrice& A) ;

    inline matrice colsref(size_t jstart, size_t jsize) const ;
    inline vecteur colref(size_t j) const;
    inline vecteur subcolref(size_t j, size_t istart, size_t isize) const ;
    inline matrice getsubmat(size_t istart, size_t isize, size_t jstart, size_t jsize) const;
    inline vecteur getcol(size_t j) const;
    inline void setcol(size_t j, const vecteur& v);
    inline vecteur getlin(size_t i) const;
    inline void setlin(size_t i, const vecteur& v);

    inline double operator[](size_t i) const ;
    inline double& operator[](size_t i) ;

    inline double operator()(size_t i,size_t j) const ;
    inline double& operator()(size_t i,size_t j) ;

    inline matrice operator()(size_t i_start, size_t i_end, size_t j_start, size_t j_end) const;

    inline const matrice& operator=(const matrice& A);
    inline const matrice& set(const double d);

    inline matrice operator*(const matrice& B) const;
    inline matrice operator*(const symmatrice& B) const;
    inline matrice operator+(const matrice& B) const;
    inline matrice operator-(const matrice& B) const;
    inline matrice operator*(double x) const;
    inline matrice operator/(double x) const;
    inline void operator+=(const matrice& B);
    inline void operator-=(const matrice& B);
    inline void operator*=(double x);
    inline void operator/=(double x);

    inline vecteur operator*(const vecteur& v) const;
    inline vecteur tmult(const vecteur &v) const;
    inline matrice tmult(const matrice &m) const;
    inline matrice multt(const matrice &m) const;
    inline matrice tmultt(const matrice &m) const;

    inline matrice transpose () const;
    inline matrice inverse() const;
    inline matrice pinverse(double reltol=0) const;
    inline void svd(matrice &U,matrice &S, matrice &V) const;
    inline double det() const;

    inline double frobenius_norm() const; 
    inline double dot(const matrice& B) const;
    inline matrice lowercholesky() const; 

    inline void write(std::ostream& f) const {
        f.write((const char*)&m,(std::streamsize)sizeof(int));
        f.write((const char*)&n,(std::streamsize)sizeof(int));
        f.write((const char*)t,(std::streamsize)(n*m*sizeof(double)));
    }

    inline void read(std::istream& f) {
        destroy();
        f.read((char*)&m,(std::streamsize)sizeof(int));
        f.read((char*)&n,(std::streamsize)sizeof(int));
        alloc(m,n);
        f.read((char*)t,(std::streamsize)(n*m*sizeof(double)));
    }

    inline void saveTxt( const char *filename ) const;
    inline void saveSubTxt( const char *filename , size_t i_start, size_t i_end, size_t j_start, size_t j_end) const;
    inline void saveBin( const char *filename ) const;
    inline void saveSubBin( const char *filename , size_t i_start, size_t i_end, size_t j_start, size_t j_end) const;
    inline void load( const char *filename );
    inline void loadTxt( const char *filename );
    inline void loadBin( const char *filename );

    friend class symmatrice;
};

inline std::ostream& operator<<(std::ostream& f,const matrice &M);

#endif

