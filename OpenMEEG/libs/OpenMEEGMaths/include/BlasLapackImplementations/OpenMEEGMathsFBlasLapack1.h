
#pragma once

//  The functions in this file are not always implemented in clapack.
//  We use their fortran versions instead,

extern "C" {
    double FC_GLOBAL(dlange,DLANGE)(const char&,const int&,const int&,const double*,const int&,double*);

    void FC_GLOBAL(dsptrf,DSPTRF)(const char&,const int&,double*,int*,int&);
    void FC_GLOBAL(dsptrs,DSPTRS)(const char&,const int&,const int&,double*,int*,double*,const int&,int&);
    void FC_GLOBAL(dsptri,DSPTRI)(const char&,const int&,double*,int*,double*,int&);
    void FC_GLOBAL(dpptrf,DPPTRF)(const char&,const int&,double*,int&);
    void FC_GLOBAL(dpptri,DPPTRI)(const char&,const int&,double*,int&);

    void FC_GLOBAL(dgesdd,DGESDD)(const char&,const int&,const int&,double*,const int&,double*,double*,const int&,double*,const int&,double*,const int&,int*,int&);
}

#define DLANGE FC_GLOBAL(dlange,DLANGE)

#define DSPTRF FC_GLOBAL(dsptrf,DSPTRF)
#define DSPTRS FC_GLOBAL(dsptrs,DSPTRS)
#define DSPTRI FC_GLOBAL(dsptri,DSPTRI)
#define DPPTRF FC_GLOBAL(dpptrf,DPPTRF)
#define DPPTRI FC_GLOBAL(dpptri,DPPTRI)

#define DGESDD FC_GLOBAL(dgesdd,DGESDD)
