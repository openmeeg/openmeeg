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

    /*
%typemap(in) const Strings* {

        std::cerr << "typemap" << std::endl;
    if(PyList_Check($input))
    {
        int size = PyList_Size($input);
        std::cerr << "typemap size " << size << std::endl;
        $1 = new Strings;
        for( int i = 0; i < size; i++)
        {
            PyObject * o = PyList_GetItem($input, i);
            std::cerr << "typemap loop " << i << std::endl;
            std::cerr << "typemap object " << o << std::endl;
            std::cerr << "typemap type is " << typeid(o).name() << std::endl;
            std::cerr << "typemap pyobject is " << Py_TYPE(o)->tp_name << std::endl;

            std::cerr << "typemap string " << PyString_AsString(o) << std::endl;

            if(PyString_Check(o)) {
                std::cerr << "typemap string " << PyString_AsString(o) << std::endl;
                $1->push_back(PyString_AsString(o));
            }
            else
            {
                PyErr_SetString(PyExc_TypeError, "list must contain strings");
                free($1);
                return NULL;
            }
        }
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "PyList of strings is expected as input");
        return NULL;
    }
}

%typemap(freearg) const Strings* {
    if ($1) {
        delete $1;
    }
}
    */

//%typemap(in) Vector& {
//     TODO check $input as np.array
//
//    PyArrayObject* input_array  = (PyArrayObject*) PyArray_FromObject($input,NPY_DOUBLE,1,1);
//
//    const size_t dim = PyArray_DIM(input_array,0);
//    std::cout << "typemap_in_Vector: " << dim << std::endl;
//
//    OpenMEEG::Vector &vec = *(new OpenMEEG::Vector(dim));
//
//    for (unsigned i=0;i<dim;++i)
//        vec(i) = *(static_cast<double*>(PyArray_GETPTR1(input_array,i)));
//
//    $1 = &vec;
//}

//%typemap(freearg) Vector& {
//    if ($1)
//        delete $1;
//}

// C++ -> Python

//%typecheck(SWIG_TYPECHECK_POINTER) Vector * {
//    $1 = ( $input != 0);
//}

//%typecheck(SWIG_TYPECHECK_VOIDPTR) Vector * {
//    $1 = ( $input != 0);
//}


//%typemap(out) Vector {
//    PyArrayObject* matarray = 0;
//
//    std::cout << "typemap_out_Vector:" << $1.size() << std::endl;
//
//    const npy_intp ndims = 1;
//    npy_intp ar_dim[] = { static_cast<npy_intp>($1.size()) };
//
//    matarray = (PyArrayObject*) PyArray_NewFromDescr(&PyArray_Type,PyArray_DescrFromType(NPY_DOUBLE),ndims,ar_dim,NULL,static_cast<void*>($1.data()),NPY_ARRAY_FARRAY,NULL);
//
//    $result = PyArray_Return(matarray);
//}

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
    Vector(PyObject* o) {
            if (!o || !PyArray_Check(o)) {
                return new Vector();
            }
            PyArrayObject *vect = (PyArrayObject *) PyArray_FromObject(o, NPY_DOUBLE, 1, 1);
            const size_t nelem = PyArray_DIM(vect, 0);

            Vector* v = new Vector(nelem);
            v->reference_data(static_cast<double *>(PyArray_GETPTR1(vect, 0)));
            return v;
    }

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
        Matrix(PyObject* o) {
            if (!o || !PyArray_Check(o)) {
                return new Matrix();
            }
            PyArrayObject *mat = (PyArrayObject *) PyArray_FromObject(o, NPY_DOUBLE, 2, 2); // able to accept a scalar, a vector and a matrix ?
            const size_t nblines = PyArray_DIM(mat, 0);
            const size_t nbcol = PyArray_DIM(mat, 1);

            Matrix* result = new Matrix(nblines, nbcol);
            result->reference_data(static_cast<double *>(PyArray_GETPTR1(mat, 0)));
            return result;
        }

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
