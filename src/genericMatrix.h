#ifndef H_genericMatrix
#define H_genericMatrix

#include <iostream>
#include "vecteur.h"

class genericMatrix
{
public:
    virtual inline int nlin() const =0;
    virtual inline int ncol() const =0;
    virtual inline double operator()(int i,int j) const =0;
    virtual inline double& operator()(int i,int j) =0;
    virtual vecteur operator*(const vecteur& x) const =0;
    virtual const genericMatrix& operator=(const double d)=0;
    virtual void saveTxt( const char *filename ) const=0;
    virtual void saveBin( const char *filename ) const=0;
    virtual void saveSubTxt( const char *filename , const int iStart, const int iEnd, const int jStart, const int jEnd) const=0;
    virtual void saveSubBin( const char *filename , const int iStart, const int iEnd, const int jStart, const int jEnd) const=0;
    virtual void loadTxt( const char *filename )=0;
    virtual void loadBin( const char *filename )=0;
    virtual void write(std::ostream& f) const =0;
    virtual void read(std::istream& f) const =0;    
};

#endif
