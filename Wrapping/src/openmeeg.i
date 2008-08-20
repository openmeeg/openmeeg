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
   #include "inversers.h"
%}

/* DLL Exports handling on Windows */
#define OPENMEEGMATHS_EXPORT
#define OPENMEEG_EXPORT

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
%include "inversers.h"

