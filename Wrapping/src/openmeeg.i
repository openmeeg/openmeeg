%module openmeeg
%{
   #include "vect3.h"
   #include "triangle.h"
   #include "vecteur.h"
   #include "vecteur.h"
   #include "matrice.h"
   #include "matrice.h"
   #include "symmatrice.h"
   #include "symmatrice.h"
   #include "sparse_matrice.h"
   #include "fast_sparse_matrice.h"
   #include "sensors.h"
   #include "geometry.h"
   #include "mesh3.h"
   #include "assemble.h"
   #include "gain.h"
   #include "forward.h"
   #include "inverse.h"
%}


%include "vect3.h"
%include "triangle.h"
%include "vecteur.h"
%include "vecteur.h"
%include "matrice.h"
%include "matrice.h"
%include "symmatrice.h"
%include "symmatrice.h"
%include "sparse_matrice.h"
%include "fast_sparse_matrice.h"
%include "geometry.h"
%include "sensors.h"
%include "mesh3.h"
%include "assemble.h"
%include "gain.h"
%include "forward.h"
%include "inverse.h"

/*%template(sparse_matrice) TsparseMatrix<double,size_t>;*/

