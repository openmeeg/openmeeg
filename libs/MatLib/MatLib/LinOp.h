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
