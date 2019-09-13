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

%inline %{

// Creator of Vector from PyArrayObject or Vector
OpenMEEG::Vector *new_OpenMEEG_Vector(PyObject *o) {
    OpenMEEG::Vector *out = nullptr;
    if (PyArray_Check(o)) {
        PyArrayObject *vect = (PyArrayObject *) PyArray_FromObject(o, NPY_DOUBLE, 1, 1);
        const size_t nelem = PyArray_DIM(vect, 0);
        OpenMEEG::Vector *v = new Vector(nelem);
        v->reference_data(static_cast<double *>(PyArray_GETPTR1(vect, 0)));
        out = v;
    } else {
        void *ptr = 0 ;
        int res = SWIG_ConvertPtr(o, &ptr, SWIGTYPE_p_OpenMEEG__Vector,  0  | 0);
        if (SWIG_IsOK(res)) {
            out = new Vector(*(reinterpret_cast<OpenMEEG::Vector *>(ptr)), DEEP_COPY);
        } else {
            PyErr_SetString(PyExc_TypeError, "Input object is neither a PyArray nor a Vector.");
            out = new OpenMEEG::Vector();
        }
    }
    return out;
}

// Creator of Vector from PyArrayObject or Matrix
OpenMEEG::Matrix *new_OpenMEEG_Matrix(PyObject *o) {
    OpenMEEG::Matrix *out = nullptr;
    if (PyArray_Check(o)) {
        int nbdims = PyArray_NDIM((PyArrayObject *)o);
        if (nbdims == 2) {
            PyArrayObject *mat = (PyArrayObject *) PyArray_FromObject(o, NPY_DOUBLE, 2, 2); // able to accept a scalar, a vector and a matrix ?
            const size_t nblines = PyArray_DIM(mat, 0);
            const size_t nbcol = PyArray_DIM(mat, 1);

            OpenMEEG::Matrix *result = new Matrix(nblines, nbcol);
            result->reference_data(static_cast<double *>(PyArray_GETPTR2(mat, 0, 0)));
            out = result;
        } else {
            PyErr_SetString(PyExc_TypeError, "Matrix requires an 2 dimensions nbarray, returning an empty matrix instead.");
            out = new Matrix();
        }
    } else {
        void *ptr = 0 ;
        int res = SWIG_ConvertPtr(o, &ptr, SWIGTYPE_p_OpenMEEG__Matrix,  SWIG_POINTER_OWN  | 0);
        if (SWIG_IsOK(res))
            out = new Matrix(*(reinterpret_cast<OpenMEEG::Matrix *>(ptr)), DEEP_COPY);
        else
            PyErr_SetString(PyExc_TypeError, "Input object is neither a PyArray nor a Matrix.");
    }
    return out;
}
%}

namespace OpenMEEG {

    // Python -> C++
%typemap(in) Vector& {
        $1 = new_OpenMEEG_Vector($input);
}

%typemap(freearg) Vector& {
        if ($1) delete $1;
}

%typemap(in) Matrix& {
        $1 = new_OpenMEEG_Matrix($input);
}

%typemap(freearg) Matrix& {
        if ($1) delete $1;
}

}

#endif // SWIGPYTHON


// /////////////////////////////////////////////////////////////////
// extensions
// /////////////////////////////////////////////////////////////////

%extend OpenMEEG::Vertex {
    // TODO almost.. if I do: v.index() I get:
    // <Swig Object of type 'unsigned int *' at 0x22129f0>
    // I want simply an unsigned. workaround:
    unsigned int getindex() {
        return ($self)->index();
    }
}

%extend OpenMEEG::Triangle{
        // TODO almost.. if I do: t.index() I get:
        // <Swig Object of type 'unsigned int *' at 0x22129f0>
        // I want simply an unsigned. workaround:
        unsigned int getindex() {
            return ($self)->index();
        }
}

%extend OpenMEEG::Vector {
    Vector(PyObject *o) {
        return new_OpenMEEG_Vector(o);
    }

    PyObject *array() {
        const npy_intp ndims = 1;
        npy_intp ar_dim[] = { static_cast<npy_intp>(($self)->size()) };
        PyArrayObject *array = (PyArrayObject*) PyArray_SimpleNewFromData(ndims, ar_dim, NPY_DOUBLE, static_cast<void*>(($self)->data()));
        return PyArray_Return(array);
    }

    // TODO almost.. v(2)=0. does not work, workaround:
    void setvalue(unsigned int i, double d) {
        (*($self))(i) = d;
    }

    double value(unsigned int i) {
        if ( i >= ($self)->size() ) {
            PyErr_SetString(PyExc_TypeError, "Out of range");
            return std::nan("");
        }
        return (*($self))(i);
    }
}

%extend OpenMEEG::Matrix {
    Matrix(PyObject *o) {
        return new_OpenMEEG_Matrix(o);
    }


    PyObject *array() {
        const npy_intp ndims = 2;
        npy_intp* ar_dim =  new npy_intp[ndims];
        ar_dim[0] = ($self)->nlin();
        ar_dim[1] = ($self)->ncol();

        PyArrayObject *array = (PyArrayObject*) PyArray_SimpleNewFromData(ndims, ar_dim, NPY_DOUBLE, static_cast<void*>(($self)->data()));
        return PyArray_Return(array);
    }

    void setvalue(unsigned int i, unsigned int j, double d) {
        (*($self))(i,j) = d;
    }

    double value(unsigned int i, unsigned int j) {
        if ( ( i >= ($self)->nlin() ) ||
             ( j >= ($self)->ncol() ) ) {
            PyErr_SetString(PyExc_TypeError, "Out of range");
            return std::nan("");
        }
        return (*($self))(i,j);
    }
}


%extend OpenMEEG::Mesh{
        // TODO almost.. if I do: m.name() I get:
        // <Swig Object of type 'std::string *' at 0x2c92ea0>
        Mesh(PyObject * py_v, PyObject * py_i) {
            if (!py_v || !PyArray_Check(py_v)) {
                return new Mesh();
            }
            PyArrayObject *mat_v = (PyArrayObject *) PyArray_FromObject(py_v, NPY_DOUBLE, 2, 2);
            int nbdims_v = PyArray_NDIM((PyArrayObject *) mat_v);
            if (nbdims_v == 2) {
                const size_t nbVertices = PyArray_DIM(mat_v, 0);
                const size_t nbcol = PyArray_DIM(mat_v, 1);
                if (nbcol < 3) {
                    PyErr_SetString(PyExc_TypeError,
                                    "Matrix of triangles requires at least 3 columns, standing for x, y, z of vertices.");
                    return new Mesh();
                } else {
                    Mesh * newMesh = new Mesh(nbVertices, nbVertices);
                    for (int cptV = 0; cptV < nbVertices; ++cptV) {
                        double x = *(double *) PyArray_GETPTR2(mat_v, cptV, 0);
                        double y = *(double *) PyArray_GETPTR2(mat_v, cptV, 1);
                        double z = *(double *) PyArray_GETPTR2(mat_v, cptV, 2);
                        Vertex * newVertex = new Vertex(x, y, z);
                        newMesh->add_vertex(*newVertex);
                    }
                    if (!py_i || !PyArray_Check(py_i)) {
                        return newMesh;
                    }
                    PyArrayObject *mat_i = (PyArrayObject *) PyArray_FromObject(py_i, NPY_DOUBLE, 2, 2);
                    int nbdims_i = PyArray_NDIM((PyArrayObject *) mat_i);
                    if (nbdims_i == 2) {
                        const size_t nbVertices = PyArray_DIM(mat_i, 0);
                        const size_t nbcol = PyArray_DIM(mat_i, 1);
                        if (nbcol != 3) {
                            PyErr_SetString(PyExc_TypeError,
                                            "Matrix of triangles requires exactly 3 columns, standing for x, y, z of indices.");
                            return newMesh;
                        } else {
                            for (int cptV = 0; cptV < nbVertices; ++cptV) {
                                int x = *(int *) PyArray_GETPTR2(mat_i, cptV, 0);
                                int y = *(int *) PyArray_GETPTR2(mat_i, cptV, 1);
                                int z = *(int *) PyArray_GETPTR2(mat_i, cptV, 2);
                                Triangle t( newMesh->vertices()[x], newMesh->vertices()[y], newMesh->vertices()[z]);
                                newMesh->push_back(t);
                            }
                            return newMesh;
                        }
                    } else {
                        PyErr_SetString(PyExc_TypeError,
                                        "Matrix of triangles requires an 2 dimensions nbarray, returning an empty matrix instead.");
                        return new Mesh();
                    }
                }
            } else {
                PyErr_SetString(PyExc_TypeError,
                                "Matrix of vertices requires an 2 dimensions nbarray, returning an empty matrix instead.");
                return new Mesh();
            }
            return new Mesh();
        }
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
