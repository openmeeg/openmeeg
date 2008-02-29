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

#ifndef H_danielsson
#define H_danielsson

#define DIM 3
#include <limits.h>
#include <float.h>
#include "vect3.h"
#include <assert.h>
#include <math.h>
#include "mesh3.h"

typedef Vect3 dpoint;
typedef Triangle ipoint;
typedef int dim_t;

const int UnknownPoint=INT_MAX;
const double UnknownDist=DBL_MAX;

double dist_point_cell(const dpoint&m ,const dpoint *pts,const ipoint& cell,dpoint& alphas,bool& inside);
double dist_point_mesh(const dpoint&m ,const Mesh &mesh,dpoint& alphas,int &nearestNumber);

#endif
