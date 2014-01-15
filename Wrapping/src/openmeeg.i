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

        /* Create a Matrix from an array
           how to assign elements, ex: A(1,2)=4. ? setvalue ? TODO */
        static OpenMEEG::Matrix fromarray(PyObject* mat) {
            if (!mat) {
                PyErr_SetString(PyExc_RuntimeError, "Zero pointer passed instead of valid array.");
                return(NULL);
            }
            PyArrayObject *matt;
            matt = (PyArrayObject *)PyArray_FromObject(mat, NPY_DOUBLE, 1, 2);
            size_t nl = matt->dimensions[0];
            size_t nc = 1;
            if (matt->nd ==2) nc = matt->dimensions[1];
            OpenMEEG::Matrix omat(nl, nc);
            for (unsigned i = 0; i< nl; ++i) {
                for (unsigned j = 0; j< nc; ++j) {
                    omat(i,j)= *(double * )(matt->data + i*matt->strides[0] + j*matt->strides[1]);
                }
            }
            return omat;
        }
        
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

namespace std {
    %template(vector_int) vector<int>;
    %template(vector_unsigned) vector<unsigned int>;
    %template(vector_double) vector<double>;
    %template(vector_triangle) vector<OpenMEEG::Triangle>;
}

namespace OpenMEEG {
    %typedef std::vector<OpenMEEG::Triangle> Triangles;
}

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

%extend OpenMEEG::Triangle {
    // TODO almost.. if I do: t.index() I get:
    // <Swig Object of type 'unsigned int *' at 0x22129f0>
    // I want simply an unsigned. workaround:
    unsigned int getindex() {
        return ($self)->index();
    }
}

%extend OpenMEEG::Vector {
    // TODO almost.. v(2)=0. does not work, workaround:
    void setvalue(unsigned int i, double d) {
        (*($self))(i)=d;
    }
}

%extend OpenMEEG::Matrix {
    void setvalue(unsigned int i, unsigned int j, double d) {
        (*($self))(i,j)=d;
    }
}
/* TODO
%include <cpointer.i>
%pointer_class(Interface ,InterfaceP)
We would like to have an Interface when asking for
i=geom.outermost_interface()
instead we have a pointer to it:
<Swig Object of type 'Interface *' at 0xa1e1590>
*/

static PyObject* asarray(OpenMEEG::Matrix* _mat);
static PyObject* asarray(OpenMEEG::Vector* _vec);
static OpenMEEG::Matrix fromarray(PyObject* _mat);
