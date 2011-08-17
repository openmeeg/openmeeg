%module(docstring="OpenMEEG bindings for python") openmeeg
%{
    #define SWIG_FILE_WITH_INIT 
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

    using namespace OpenMEEG;

    #ifdef SWIGPYTHON

        #include <numpy/arrayobject.h>

        static PyObject* asarray(Matrix* _mat) {
            if (!_mat) {
                PyErr_SetString(PyExc_RuntimeError, "Zero pointer passed instead of valid Matrix struct.");
                return(NULL);
            }
            
            /* array object */
            PyArrayObject* matarray = 0;
            
            /* Get the number of dimensions from the Matrix
             */
            int ar_dim[7], ndims;
            ndims = 2;
            ar_dim[0] = _mat->ncol();
            ar_dim[1] = _mat->nlin();
            
            /* create numpy array */
            matarray = (PyArrayObject*) PyArray_FromDimsAndData(ndims,ar_dim,PyArray_DOUBLE,(char*) _mat->data());
            int tmp = matarray->strides[0];
            matarray->strides[0] = matarray->strides[1];
            matarray->strides[1] = tmp;
            
            tmp = matarray->dimensions[0];
            matarray->dimensions[0] = matarray->dimensions[1];
            matarray->dimensions[1] = tmp;
            
            return PyArray_Return((PyArrayObject*) matarray);
        }

        static PyObject* asarray(Vector* _vec) {
            if (!_vec) {
                PyErr_SetString(PyExc_RuntimeError, "Zero pointer passed instead of valid Vector struct.");
                return(NULL);
            }
            
            /* array object */
            PyArrayObject* matarray = 0;
            
            /* Get the size of the Vector
             */
            int ar_dim[7], ndims;
            ndims = 1;
            ar_dim[0] = _vec->size();
            
            /* create numpy array */
            matarray = (PyArrayObject*) PyArray_FromDimsAndData (ndims,ar_dim,PyArray_DOUBLE,(char*) _vec->data());
            
            return PyArray_Return ((PyArrayObject*) matarray);
        }

    #endif
%}

%include "numpy.i"

%init %{ 
import_array(); 
%}

%exception {
    try {
        $action
    }
    catch (std::exception& e) {
        PyErr_SetString(PyExc_RuntimeError,e.what());
        return NULL;
    }
}

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

static PyObject* asarray(Matrix* _mat);
static PyObject* asarray(Vector* _vec);
