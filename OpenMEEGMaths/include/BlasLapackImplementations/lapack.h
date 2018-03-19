#pragma once

extern "C" {
    void LAPACK(dgetrf,DGETRF)(const int&,const int&,double*,const int&,int*,int&);
    void LAPACK(dgetri,DGETRI)(const int&,double*,const int&,int*,double*,const int&,int&);
}
