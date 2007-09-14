#ifndef H_LINOP
#define H_LINOP

#include "vecteur_dcl.h"

class LinOp
{
public:
    virtual vecteur operator * ( const vecteur &x) const = 0;
    virtual ~LinOp() {};
};

#endif
