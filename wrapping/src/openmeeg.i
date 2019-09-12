%module(docstring="OpenMEEG bindings for python") openmeeg

%include <exception.i>
%exception {
    try {
        $action
    } catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
}

#ifdef SWIGWIN
%include <windows.i>
#endif



%include <std_string.i>
%include <std_vector.i>

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
    #include <interface.h>
    #include <domain.h>
    #include <assemble.h>
    #include <gain.h>
    #include <forward.h>
    #include <iostream>

    #ifdef SWIGPYTHON

        #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
        #include <numpy/arrayobject.h>

    #endif

    using namespace OpenMEEG;
%}

// /////////////////////////////////////////////////////////////////
// Preprocessing setup
// /////////////////////////////////////////////////////////////////

#pragma SWIG nowarn=302, 315, 389, 401, 509, 801, 472, 473, 476, 362, 503, 514, 516, 842, 845

// /////////////////////////////////////////////////////////////////
// Ignore rules for operators
// /////////////////////////////////////////////////////////////////

%ignore operator>>;
%ignore operator<<;
%ignore operator==;
%ignore operator[];
%ignore operator!=;
%ignore operator*=;
%ignore operator/=;
//%ignore operator bool;
//%ignore operator int;
//%ignore operator float;
//%ignore operator double;
//%ignore operator double *;

// /////////////////////////////////////////////////////////////////
// Legacy
// /////////////////////////////////////////////////////////////////

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

// /////////////////////////////////////////////////////////////////
// Definitions
// /////////////////////////////////////////////////////////////////

#define OPENMEEGMATHS_EXPORT
#define OPENMEEG_EXPORT


namespace std {
    %template(vector_int) vector<int>;
    %template(vector_unsigned) vector<unsigned int>;
    %template(vector_double) vector<double>;
    %template(vector_vertex) vector<OpenMEEG::Vertex>;
    %template(vector_pvertex) vector<OpenMEEG::Vertex *>;
    %template(vector_triangle) vector<OpenMEEG::Triangle>;
    %template(vector_mesh) vector<OpenMEEG::Mesh>;
    %template(vector_string) vector<std::string>;
    %template(vector_interface) vector<OpenMEEG::Interface>;
}

namespace OpenMEEG {
    %typedef std::vector<OpenMEEG::Vertex> Vertices;
    %typedef std::vector<OpenMEEG::Vertex *> PVertices;
    %typedef std::vector<OpenMEEG::Triangle> Triangles;
    %typedef std::vector<OpenMEEG::Mesh> Meshes;
    %typedef std::vector<std::string>    Strings;
}

// /////////////////////////////////////////////////////////////////
// Typemaps
// /////////////////////////////////////////////////////////////////

#ifdef SWIGPYTHON

namespace OpenMEEG {


    // Python -> C++
%typemap(in) Vector& {
        if (PyArray_Check($input)) {
            PyArrayObject *vect = (PyArrayObject *) PyArray_FromObject($input, NPY_DOUBLE, 1, 1);
            const size_t nelem = PyArray_DIM(vect, 0);
            Vector *v = new Vector(nelem);
            v->reference_data(static_cast<double *>(PyArray_GETPTR1(vect, 0)));
            $1 = v;
        } else {
            void *ptr = 0 ;
            int res = SWIG_ConvertPtr($input, &ptr, SWIGTYPE_p_OpenMEEG__Vector,  0  | 0);
            if (SWIG_IsOK(res))
                $1 = reinterpret_cast<OpenMEEG::Vector *>(ptr);
            else
                PyErr_SetString(PyExc_TypeError, "Input object is neither a PyArray nor a Vector.");
        }
}

%typemap(in) Matrix& {
            int nbdims = PyArray_NDIM((PyArrayObject *)$input);
            if (nbdims == 2) {
                PyArrayObject *mat = (PyArrayObject *) PyArray_FromObject($input, NPY_DOUBLE, 2,
                                                                          2); // able to accept a scalar, a vector and a matrix ?
                const size_t nblines = PyArray_DIM(mat, 0);
                const size_t nbcol = PyArray_DIM(mat, 1);

                Matrix * result = new Matrix(nblines, nbcol);
                result->reference_data(static_cast<double *>(PyArray_GETPTR2(mat, 0, 0)));
                $1 = result;
            } else {
                PyErr_SetString(PyExc_TypeError, "Matrix requires an 2 dimensions nbarray, returning an empty matrix instead.");
                $1 = new Matrix();
            }
}

}

#endif // SWIGPYTHON

// /////////////////////////////////////////////////////////////////
// extensions
// /////////////////////////////////////////////////////////////////

%extend OpenMEEG::Sensors {
    Sensors(const char*, Vector& a) {
        Sensors *S = new Sensors();
        //
        std::cout << "extend_sensors: " << a << std::endl;
        return S;
    }
}

%extend OpenMEEG::Vertex {
    // TODO almost.. if I do: v.index() I get:
    // <Swig Object of type 'unsigned int *' at 0x22129f0>
    // I want simply an unsigned. workaround:
    unsigned int getindex() {
        return ($self)->index();
    }
}

%extend OpenMEEG::Triangle {
    // TODO almost.. if I do: t.index() I get:
    // <Swig Object of type 'unsigned int *' at 0x22129f0>
    // I want simply an unsigned. workaround:
    unsigned int getindex() {
        return ($self)->index();
    }
}

%extend OpenMEEG::Vector {

    PyObject* array() {
            const npy_intp ndims = 1;
            npy_intp ar_dim[] = { static_cast<npy_intp>(($self)->size()) };

            PyArrayObject* array = (PyArrayObject*) PyArray_SimpleNewFromData(ndims, ar_dim, NPY_DOUBLE, static_cast<void*>(($self)->data()));
            return PyArray_Return(array);
    }

    // TODO almost.. v(2)=0. does not work, workaround:
    void setvalue(unsigned int i, double d) {
        (*($self))(i)=d;
    }
}

%extend OpenMEEG::Matrix {
        PyObject* array() {
            const npy_intp ndims = 2;
            npy_intp* ar_dim =  new npy_intp[ndims];
            ar_dim[0] = ($self)->nlin();
            ar_dim[1] = ($self)->ncol();

            PyArrayObject* array = (PyArrayObject*) PyArray_SimpleNewFromData(ndims, ar_dim, NPY_DOUBLE, static_cast<void*>(($self)->data()));
            return PyArray_Return(array);
        }

        void setvalue(unsigned int i, unsigned int j, double d) {
        (*($self))(i,j)=d;
    }
}


%extend OpenMEEG::Mesh {
    // TODO almost.. if I do: m.name() I get:
    // <Swig Object of type 'std::string *' at 0x2c92ea0>
    const char* __str__() {
        return ($self)->name().c_str();
    }
}

// /////////////////////////////////////////////////////////////////
// Input
// /////////////////////////////////////////////////////////////////

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
%include <interface.h>
%include <domain.h>
%include <assemble.h>
%include <gain.h>
%include <forward.h>

// /////////////////////////////////////////////////////////////////
//
// /////////////////////////////////////////////////////////////////
