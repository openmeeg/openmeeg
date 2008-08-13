%module(docstring="OpenMEEG bidings for python") openmeeg
%{
   #include "vect3.h"
   #include "triangle.h"
   #include "linop.h"
   #include "vector.h"
   #include "matrix.h"
   #include "symmatrix.h"
   #include "sparse_matrix.h"
   #include "fast_sparse_matrix.h"
   #include "sensors.h"
   #include "geometry.h"
   #include "mesh3.h"
   #include "assemble.h"
   #include "gain.h"
   #include "forward.h"
   #include "inverse.h"
%}

#ifdef DOCSTRINGS
%include "docstrings.i"
#endif

%include "vect3.h"
%include "triangle.h"
%include "linop.h"
%include "vector.h"
%include "matrix.h"
%include "symmatrix.h"
%include "sparse_matrix.h"
%include "fast_sparse_matrix.h"
%include "geometry.h"
%include "sensors.h"
%include "mesh3.h"
%include "assemble.h"
%include "gain.h"
%include "forward.h"
%include "inverse.h"

