#ifndef H_Linop
#define H_Linop

#include "vecteur_dcl.h"

class LinOp
{
public:
    virtual vecteur operator * ( const vecteur &x) const = 0;
    virtual ~LinOp() {};
};

#endif
