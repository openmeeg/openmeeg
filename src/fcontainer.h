#ifndef H_fContainer
#define H_fContainer

#include <math.h>
#include "vect3.h"

template<class T> class fContainer
{
    // fContainer class , an interface for a class having a T-valued function in \real^3
public:
    fContainer(){}
    virtual ~fContainer(){}

    virtual T f(const Vect3& v) const=0;
};
#endif

