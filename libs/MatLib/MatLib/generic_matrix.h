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

#ifndef H_genericMatrix

#define H_genericMatrix
#include <iostream>
#include "LinOp.h"

template<class T> class TgenericMatrix : public LinOp
{
public:
    virtual size_t nlin() const=0;
    virtual size_t ncol() const=0;
    virtual T operator()(size_t i,size_t j) const=0;
    virtual T& operator()(size_t i,size_t j)=0;
    virtual void saveTxt( const char *filename ) const=0;
    virtual void saveBin( const char *filename ) const=0;
    virtual void loadTxt( const char *filename )=0;
    virtual void loadBin( const char *filename )=0;
    virtual void write(std::ostream& f) const =0;
    virtual void read(std::istream& f) =0;
    virtual std::ostream& operator>>(std::ostream& f) const =0;

// #ifdef USE_MATIO
//     virtual void saveMat( const char *filename ) const=0;
//     virtual void loadMat( const char *filename )=0;
// #endif
    
};

typedef TgenericMatrix<double> genericMatrix;


#endif

