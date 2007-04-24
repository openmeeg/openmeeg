#ifndef H_danielsson
#define H_danielsson

#define DIM 3
#include <limits.h>
#include <float.h>
#include "vect3.h"
#include <assert.h>
#include <math.h>
#include "mesh3.h"
typedef vect3 dpoint;
typedef triangle ipoint;
typedef int dim_t;

const int UnknownPoint=INT_MAX;
const double UnknownDist=DBL_MAX;



double dist_point_cell(const dpoint&m ,const dpoint *pts,const ipoint& cell,dpoint& alphas,bool& inside);
double dist_point_mesh(const dpoint&m ,const mesh &mesh,dpoint& alphas,int &nearestNumber);

#endif
