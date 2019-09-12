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
}

// /////////////////////////////////////////////////////////////////
// Typemaps
// /////////////////////////////////////////////////////////////////

#ifdef SWIGPYTHON

namespace OpenMEEG {


// Python -> C++

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
