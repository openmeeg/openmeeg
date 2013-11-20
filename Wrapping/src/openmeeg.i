%module(docstring="OpenMEEG bindings for python") openmeeg

%include "std_string.i"
%include "std_vector.i"

%{
    #define SWIG_FILE_WITH_INIT
    #include <vect3.h>
    #include <vertex.h>
    #include <triangle.h>
    #include <linop.h>
    #include <vector.h>
    #include <matrix.h>
    #include <symmatrix.h>
    #include <sparse_matrix.h>
    #include <fast_sparse_matrix.h>
    #include <sensors.h>
    #include <geometry.h>
    #include <geometry_io.h>
    #include <mesh.h>
    #include <domain.h>
    #include <interface.h>
    #include <assemble.h>
    #include <gain.h>
    #include <forward.h>

    using namespace OpenMEEG;

    #ifdef SWIGPYTHON

        #include <numpy/arrayobject.h>

        static PyObject* asarray(OpenMEEG::Matrix* _mat) {
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

        static PyObject* asarray(OpenMEEG::Vector* _vec) {
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

        /* Create a Matrix or a Vector from an array
           how to assign elements, ex: A(1,2)=4. ?
        TODO:   check http://docs.scipy.org/doc/numpy/reference/c-api.types-and-structures.html at PyArray_GetPtr
        static OpenMEEG::Vector* fomarray(PyArrayObject* vec) {
            if (!vec) {
                PyErr_SetString(PyExc_RuntimeError, "Zero pointer passed instead of valid array.");
                return(NULL);
            }
            // create OM vector 
            OpenMEEG::Vector vect(PyArray_Size(vec));
            vect.data() = PyArray_GetPtr(vec,1); // this is not the correct way
            return vect;
        }
        */

    #endif
%}

%pythoncode {
def loadmat(fname):
    try:
        from scipy import io
        io.loadmat(fname)['linop']
    except:
        import h5py
        return h5py.File(fname)['linop'].value.T
}

%include "numpy.i"

%init %{
import_array();
//import_array1(NULL); For future python 3.x
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
%include <docstrings.i>
#endif

%include <vect3.h>
%include <vertex.h>
%include <triangle.h>
%include <linop.h>
%include <vector.h>
%include <matrix.h>
%include <symmatrix.h>
%include <sparse_matrix.h>
%include <fast_sparse_matrix.h>
%include <geometry.h>
%include <geometry_io.h>
%include <sensors.h>
%include <mesh.h>
%include <domain.h>
%include <interface.h>
%include <assemble.h>
%include <gain.h>
%include <forward.h>

/* TODO
%include <cpointer.i>
%pointer_class(Interface ,InterfaceP)
We would like to have an Interface when asking for
i=geom.outermost_interface()
instead we have:
<Swig Object of type 'Interface *' at 0xa1e1590>
*/

static PyObject* asarray(OpenMEEG::Matrix* _mat);
static PyObject* asarray(OpenMEEG::Vector* _vec);
