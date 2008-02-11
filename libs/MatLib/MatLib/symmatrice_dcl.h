#ifndef H_SYMMATRICE_DCL
#define H_SYMMATRICE_DCL

#include "generic_matrix.h"
#include <fstream>
#include <cassert>

#include "vecteur_dcl.h"

class matrice;
class vecteur;
class symmatrice;

class symmatrice : public genericMatrix  {
    size_t n;
    double *t;
    int *count;

    inline void alloc(size_t N);
    inline void copy(const symmatrice& A);
    inline void destroy();
public:
    inline symmatrice();
    inline symmatrice(size_t N) ;
    inline symmatrice(const char *filename, const char c='t');
    inline symmatrice(const symmatrice& A);
    inline explicit symmatrice(const vecteur& v);
    inline explicit symmatrice(const matrice& A); // upper triangle based (lower triangle need not to be set)
    inline symmatrice(double* T, int* COUNT, size_t N);
    inline  ~symmatrice() { destroy(); }
    inline symmatrice duplicate() const;
    inline size_t nlin() const ;
    inline size_t ncol() const ;
    inline bool empty() const ;
    inline void set(double x) ;
    inline double* DangerousGetData () ;
    inline int* DangerousGetCount () ;

    inline double operator()(size_t i,size_t j) const;
    inline double& operator()(size_t i,size_t j) ;

    inline matrice operator()(size_t i_start, size_t i_end, size_t j_start, size_t j_end) const;
    inline matrice getsubmat(size_t istart, size_t isize, size_t jstart, size_t jsize) const;
    inline symmatrice getsubmat(size_t istart, size_t iend) const;
    inline vecteur solveLin(const vecteur &B) const;
    inline void solveLin(vecteur * B, int nbvect);

    inline const symmatrice& operator=(const double d);
    inline const symmatrice& operator=(const symmatrice& A);

    inline symmatrice operator+(const symmatrice& B) const;
    inline symmatrice operator-(const symmatrice& B) const;
    inline symmatrice operator*(double x) const;
    inline symmatrice operator/(double x) const {return (*this)*(1/x);}
    inline void operator +=(const symmatrice& B);
    inline void operator -=(const symmatrice& B);
    inline void operator *=(double x);
    inline void operator /=(double x) ;
    inline matrice operator*(const matrice& B) const; // faux !!
    inline vecteur operator*(const vecteur& v) const; // faux ?

    inline symmatrice inverse() const;
    inline symmatrice posdefinverse() const;
    inline double det();
    inline void eigen(matrice & Z, vecteur & D );

    inline void write(std::ostream& f) const ;

    inline void read(std::istream& f) ;

    inline void saveSubTxt( const char *filename , size_t i_start, size_t i_end, size_t j_start, size_t j_end) const;
    inline void saveSubBin( const char *filename , size_t i_start, size_t i_end, size_t j_start, size_t j_end) const;
    inline void saveTxt( const char *filename ) const;
    inline void saveBin( const char *filename ) const;
    inline void save( const char *filename ) const;
    inline void load( const char *filename );
    inline void loadTxt( const char *filename );
    inline void loadBin( const char *filename );
    
    inline std::ostream& operator>>(std::ostream& f) const {
        for (size_t i=0;i<this->nlin();i++) {
            for (size_t j=0;j<this->ncol();j++) {
                f << (*this)(i,j) << " ";
            }
            f << std::endl;
        }
        return f;
    }

    friend class matrice;
};

#endif

