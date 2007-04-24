#ifndef H_fast_sparse_matrice_dcl
#define H_fast_sparse_matrice_dcl

#include "sparse_matrice.h"

class vecteur;

class fast_sparse_matrice : public genericMatrix
{
public:
    typedef sparse_matrice::idxType idxType;
    typedef sparse_matrice::valType valType;

protected:

    valType *tank;
    idxType *js;
    idxType *rowindex;
    idxType nlignes;
    idxType ncolonnes;

    inline void alloc(idxType nl, idxType nc, idxType nz);
    inline void destroy();

public:
    inline fast_sparse_matrice();
    inline fast_sparse_matrice( const sparse_matrice &M);
    inline fast_sparse_matrice( const fast_sparse_matrice &M);
    inline ~fast_sparse_matrice() {destroy();}
    inline size_t nlin() const ;
    inline size_t ncol() const ;
    inline void saveTxt( const char *filename ) const;
    inline void saveBin( const char *filename ) const;
    inline void loadTxt( const char *filename );
    inline void loadBin( const char *filename );
    inline void write(std::ostream& f) const;
    inline void read(std::istream& f);

    inline double* getbuf() const;
    inline double operator()(size_t i,size_t j) const;
    inline double& operator()(size_t i,size_t j);
    inline vecteur operator * (const vecteur &v) const;
    inline void operator =( const fast_sparse_matrice &M);
    inline vecteur solve(const vecteur &rhs_vec) const;

    inline friend std::ostream& operator<<(std::ostream& f,const fast_sparse_matrice &M);
};

#endif
