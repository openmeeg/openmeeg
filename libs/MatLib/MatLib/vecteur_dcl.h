#ifndef H_vecteur_dcl
#define H_vecteur_dcl

#include "MatLibConfig.h"
#include <iostream>
#include <cassert>

class matrice;
class symmatrice;

class vecteur {
    size_t n;
    double *t;
    int *count;

    inline void alloc(size_t N);
    inline void destroy();
    inline void copy(const vecteur& A) ;
    inline void copyout(double * p) const;
    inline void copyin(const double * p);
public:
    inline vecteur();
    inline vecteur(size_t N);
    inline vecteur(const vecteur& A);
    inline explicit vecteur(matrice& A);
    inline explicit vecteur(symmatrice& A);
    inline vecteur(double* T, int* COUNT, size_t N);
    inline  ~vecteur() { destroy(); }

    inline size_t size() const ;
    inline bool empty() const ;
    inline void DangerousBuild(double *t, size_t size);
    inline void DangerousKill();
    inline double* DangerousGetData () const ;
    inline int* DangerousGetCount () const ;

    inline vecteur duplicate() const ;
    inline void copyin(const vecteur& A) ;

    inline vecteur subvectref(size_t istart, size_t size) const ;
    inline double operator()(size_t i) const ;
    inline double& operator()(size_t i) ;
    inline const vecteur& operator=(const vecteur& A);

    inline vecteur operator+(const vecteur& v) const;
    inline vecteur operator-(const vecteur& v) const;
    inline void operator+=(const vecteur& v);
    inline void operator-=(const vecteur& v);
    inline void operator*=(double x);
    inline void operator/=(double x);
    inline vecteur operator+(double i) const;
    inline vecteur operator-(double i) const;
    inline vecteur operator*(double x) const;
    inline vecteur operator/(double x) const ;
    inline double operator*(const vecteur& v) const;

    inline vecteur conv(const vecteur& v) const;
    inline vecteur conv_trunc(const vecteur& v) const;
    inline matrice outer_product(const vecteur& v) const;

    inline double norm()const;
    inline size_t argabsmin() const;
    inline size_t argabsmax() const;
    inline double abssum() const;
    inline double sum() const;
    inline double mean() const;

    inline vecteur div(const vecteur &d) const;
    inline vecteur invsqrt() const;
    inline vecteur sqrt() const;
    inline vecteur inv() const;
    inline vecteur sin() const;
    inline vecteur cos() const;
    inline vecteur exp() const;
    inline vecteur ln() const;
    inline vecteur log10() const;
    inline vecteur pow( double x) const;
    inline vecteur pow( const vecteur &x) const;
    inline vecteur erf() const;
    inline vecteur erfc() const;

    inline void set(double x);
    inline void write(std::ostream& f) const ;
    inline void read(std::istream& f) ;
    inline void saveTxt( const char* filename) const;
    inline void saveBin( const char* filename) const;
    inline void loadTxt( const char* filename);
    inline void loadBin( const char* filename);

    friend class symmatrice;
    friend class matrice;
};

inline vecteur operator * (const double &d, const vecteur &v) ;
inline std::ostream& operator<<(std::ostream& f,const vecteur &M);
inline std::istream& operator>>(std::istream& f,vecteur &M);

#endif

