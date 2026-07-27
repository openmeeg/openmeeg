#pragma once

#include <BlasLapackImplementations/FortranAlias.h>

extern "C" {
    void LAPACK(dgetrf,DGETRF)(const int&,const int&,double*,const int&,int*,int&) FC_ALIAS(dgetrf);
    void LAPACK(dgetri,DGETRI)(const int&,double*,const int&,int*,double*,const int&,int&) FC_ALIAS(dgetri);
}
