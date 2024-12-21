%module(docstring="OpenMEEG bindings for python") openmeeg

// TODO: Should use modern macros https://numpy.org/doc/stable/reference/swig.interface-file.html#macros

%feature("autodoc", "1");

%inline %{

    #include <string>
    #include <exception>

    class Error: public std::exception {
    public:

        Error(const int t,const char* m): typ(t),mes(m) { }

        int type() const { return typ; }

        const std::string& message() const { return mes; }

    private:

        const int         typ;
        const std::string mes;
    };
%}

%include <exception.i>
%exception {
    try {
        $action
    } catch (const Error& e) {
        SWIG_exception(e.type(),e.message().c_str());
    } catch (const OpenMEEG::maths::IOException & e) {
        SWIG_exception(SWIG_IOError,e.what());
    } catch (const OpenMEEG::IOException & e) {
        SWIG_exception(SWIG_IOError,e.what());
    } catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError,e.what());
    }
}

#define SWIG_PYTHON_SILENT_MEMLEAK

#ifdef SWIGWIN
%include <windows.i>
#endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4127)  /* conditional expression is constant */
#pragma warning(disable:4701)  /* potentially uninitialized local variable 'v' used */
#endif

%include <std_string.i>
%include <std_vector.i>

%{
    #include <logger.h>
    #include <vect3.h>
    #include <vertex.h>
    #include <triangle.h>
    #include <linop.h>
    #include <vector.h>
    #include <matrix.h>
    #include <symmatrix.h>
    #include <sparse_matrix.h>
    #include <block_matrix.h>
    #include <fast_sparse_matrix.h>
    #include <sensors.h>
    #include <geometry.h>
    #include <GeometryIO.h>
    #include <mesh.h>
    #include <integrator.h>
    #include <interface.h>
    #include <domain.h>
    #include <assemble.h>
    #include <gain.h>
    #include <forward.h>
    #include <iostream>

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

%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%{
#undef SWIG_FILE_WITH_INIT
%}

%init %{
    import_array();
%}

// /////////////////////////////////////////////////////////////////
// Definitions
// /////////////////////////////////////////////////////////////////

#define OPENMEEGMATHS_EXPORT
#define OPENMEEG_EXPORT

namespace std {
    // std::vector<Mesh> cannot be handled similarly because swig assumes a copy constructor.
    %template(vector_int) vector<int>;
    %template(vector_unsigned) vector<unsigned int>;
    %template(vector_double) vector<double>;
    %template(vector_vertex) vector<OpenMEEG::Vertex>;
    %template(vector_pvertex) vector<OpenMEEG::Vertex *>;
    %template(vector_triangle) vector<OpenMEEG::Triangle>;
    %template(vector_string) vector<std::string>;
    %template(vector_interface) vector<OpenMEEG::Interface>;
    %template(vector_simple_dom) vector<OpenMEEG::SimpleDomain>;
    %template(vector_domain) vector<OpenMEEG::Domain>;
    %template(vector_oriented_mesh) vector<OpenMEEG::OrientedMesh>;
}

namespace OpenMEEG {
    // The OpenMEEG:: prefix seems required...
    // Otherwise some wrapping tests are failing.
    %typedef std::vector<OpenMEEG::Vertex>   Vertices;
    %typedef std::vector<OpenMEEG::Vertex*>  PVertices;
    %typedef std::vector<OpenMEEG::Triangle> Triangles;
    %typedef std::vector<OpenMEEG::Domain>   Domains;
    %typedef std::vector<std::string>        Strings;
    %typedef std::vector<SimpleDomain>       Boundaries;
    %typedef std::vector<OrientedMesh>       OrientedMeshes;
}

namespace OpenMEEG {

    %naturalvar Logger;
    class Logger;

    %naturalvar Mesh;
    class Mesh;

    %naturalvar Matrix;
    class Matrix;

    %naturalvar SymMatrix;
    class SymMatrix;

    %naturalvar Triangle;
    class Triangle;

    %naturalvar Vect3;
    class Vect3;

    %naturalvar Vertex;
    class Vertex;

    %naturalvar Vertices;
    class Vertices;

    %naturalvar Vector;
    class Vector;

    %naturalvar Domains;
    class Domains;

    %naturalvar Boundaries;
    class Boundaries;

    %naturalvar OrientedMeshes;
    class OrientedMeshes;
}

%inline %{

    typedef std::shared_ptr<double[]> SharedData;

    // Creator of Vector from PyArrayObject or Vector

    OpenMEEG::Vector* new_OpenMEEG_Vector(PyObject* pyobj) {
        if (pyobj && PyArray_Check(pyobj)) {
            PyArrayObject* vect = reinterpret_cast<PyArrayObject*>(PyArray_FromObject(pyobj,NPY_DOUBLE,1,1));
            const size_t nelem = PyArray_DIM(vect,0);
            OpenMEEG::Vector* v = new Vector(nelem);
            v->reference_data(static_cast<double*>(PyArray_GETPTR1(vect,0)));
            return v;
        }

        //  If the object is an OpenMEEG vector converted to python, copy the vector.
        //  TODO: do we need this ???

        void* ptr = 0 ;
        if (!SWIG_IsOK(SWIG_ConvertPtr(pyobj,&ptr,SWIGTYPE_p_OpenMEEG__Vector,SWIG_POINTER_EXCEPTION)))
            throw Error(SWIG_TypeError, "Input object is neither a PyArray nor a Vector.");

        return new Vector(*(reinterpret_cast<OpenMEEG::Vector*>(ptr)),DEEP_COPY);
    }

    // Creator of Matrix from PyArrayObject or Matrix

    OpenMEEG::Matrix* new_OpenMEEG_Matrix(PyObject* pyobj) {
        if (pyobj && PyArray_Check(pyobj)) {
            const int nbdims = PyArray_NDIM(reinterpret_cast<PyArrayObject*>(pyobj));
            if (nbdims!=2)
                throw Error(SWIG_TypeError, "Matrix can only have 2 dimensions.");

            PyArrayObject* mat = reinterpret_cast<PyArrayObject*>(PyArray_FromObject(pyobj,NPY_DOUBLE,2,2));

            if (!PyArray_ISFARRAY(mat))
                throw Error(SWIG_TypeError, "Matrix requires the use of Fortran order.");

            const size_t nblines = PyArray_DIM(mat,0);
            const size_t nbcol   = PyArray_DIM(mat,1);

            OpenMEEG::Matrix* result = new Matrix(nblines,nbcol);
            result->reference_data(static_cast<double*>(PyArray_GETPTR2(mat,0,0)));
            return result;
        }

        //  If the object is an OpenMEEG matrix converted to python, just return the matrix.

        void* ptr = 0;
        if (!SWIG_IsOK(SWIG_ConvertPtr(pyobj,&ptr,SWIGTYPE_p_OpenMEEG__Matrix,SWIG_POINTER_EXCEPTION)))
            throw Error(SWIG_TypeError, "Input object must be a PyArray or an OpenMEEG Matrix.");

        return new Matrix(*(reinterpret_cast<OpenMEEG::Matrix*>(ptr)));
    }

    OpenMEEG::SymMatrix* new_OpenMEEG_SymMatrix(PyObject* pyobj) {
        if (pyobj && PyArray_Check(pyobj)) {
            const int nbdims = PyArray_NDIM(reinterpret_cast<PyArrayObject*>(pyobj));
            if (nbdims!=1)
                throw Error(SWIG_TypeError, "SymMatrix are stored as 1 dimensional arrays.");

            PyArrayObject* mat = reinterpret_cast<PyArrayObject*>(PyArray_FromObject(pyobj,NPY_DOUBLE,1,1));

            if (!PyArray_ISFARRAY(mat))
                throw Error(SWIG_TypeError, "SymMatrixMatrix requires the use of Fortran order.");

            const size_t size = PyArray_DIM(mat,0);
            const size_t nblines = round((-1 + sqrt(1+8*size))/2);

            OpenMEEG::SymMatrix* result = new SymMatrix(nblines);
            result->reference_data(static_cast<double*>(PyArray_GETPTR2(mat,0,0)));
            return result;
        }

        //  If the object is an OpenMEEG matrix converted to python, just return the matrix.

        void* ptr = 0;
        if (!SWIG_IsOK(SWIG_ConvertPtr(pyobj,&ptr,SWIGTYPE_p_OpenMEEG__Matrix,SWIG_POINTER_EXCEPTION)))
            throw Error(SWIG_TypeError, "Input object must be a PyArray or an OpenMEEG Matrix.");

        return new SymMatrix(*(reinterpret_cast<OpenMEEG::SymMatrix*>(ptr)));
    }

    IndexMap
    geom_add_vertices(Geometry* geom,PyObject* pyobj) {

        if (pyobj==nullptr || !PyArray_Check(pyobj))
            throw Error(SWIG_TypeError,"Vertices matrix should be an array.");

        PyArrayObject* array = reinterpret_cast<PyArrayObject*>(PyArray_FromObject(pyobj,NPY_DOUBLE,0,0));
        if (array==nullptr)
            throw Error(SWIG_ValueError,"Vertices matrix cannot be converted into a matrix of double.");

        if (PyArray_NDIM(array)!=2 || PyArray_DIM(array,1)!=3)
            throw Error(SWIG_ValueError,"Vertices matrix must be a 2 dimensional array with 3 columns.");

        IndexMap indmap;
        const size_t num_vertices = PyArray_DIM(array,0);
        for (unsigned int i=0; i<num_vertices; ++i) {
            const double x = *reinterpret_cast<double*>(PyArray_GETPTR2(array,i,0));
            const double y = *reinterpret_cast<double*>(PyArray_GETPTR2(array,i,1));
            const double z = *reinterpret_cast<double*>(PyArray_GETPTR2(array,i,2));
            indmap.insert({ i, geom->add_vertex(Vertex(x,y,z)) });
        }
        return indmap;
    }

    void
    mesh_add_triangles(Mesh* mesh,PyObject* pyobj,const IndexMap& indmap) {
        if (pyobj==nullptr || !PyArray_Check(pyobj))
            throw Error(SWIG_TypeError,"Matrix of triangles should be an array.");

        PyArrayObject* array = reinterpret_cast<PyArrayObject*>(pyobj);
        if (PyArray_SIZE(array)==0) {
            std::ostringstream oss;
            oss << "Matrix of triangles for mesh \"" << mesh->name() << "\" was empty";
            throw Error(SWIG_ValueError,oss.str().c_str());
        }
        const PyArray_Descr* descr = PyArray_DESCR(array);
        const int type_num = descr->type_num;
        if (!PyArray_EquivTypenums(type_num,NPY_INT32) &&
            !PyArray_EquivTypenums(type_num,NPY_UINT32) &&
            !PyArray_EquivTypenums(type_num,NPY_INT64) &&
            !PyArray_EquivTypenums(type_num,NPY_UINT64)) {
            std::ostringstream oss;
            oss << "Wrong dtype for triangles array (only 32 or 64 int or uint supported), got type '" << descr->kind << PyDataType_ELSIZE(descr) << "'";
            throw Error(SWIG_TypeError, oss.str().c_str());
        }

        const size_t ndims = PyArray_NDIM(array);
        if (ndims!=2)
            throw Error(SWIG_TypeError,"Matrix of triangles must be a 2 dimensional array.");

        const size_t nbTriangles  = PyArray_DIM(array,0);
        if (PyArray_DIM(array,1)!=3)
            throw Error(SWIG_TypeError,"Matrix of triangles requires exactly 3 columns, standing for indices of 3 vertices.");

        mesh->reference_vertices(indmap);

        auto get_vertex = [&](PyArrayObject* mat,const int i,const int j) {
            const unsigned vi = *reinterpret_cast<unsigned*>(PyArray_GETPTR2(mat,i,j));
            if (vi>=indmap.size()) {
                std::ostringstream oss;
                oss << "Vertex index " << vi << " in triangle " << i << " out of range";
                throw Error(SWIG_ValueError,oss.str().c_str());
            }
            return &(mesh->geometry().vertices().at(indmap.at(vi)));
        };

        for (int unsigned i=0; i<nbTriangles; ++i) {
            Vertex* v1 = get_vertex(array,i,0);
            Vertex* v2 = get_vertex(array,i,1);
            Vertex* v3 = get_vertex(array,i,2);
            mesh->triangles().push_back(Triangle(v1,v2,v3));
        }
    }
%}

// /////////////////////////////////////////////////////////////////
// Typemaps
// /////////////////////////////////////////////////////////////////

namespace OpenMEEG {

    // Python -> C++
    %typecheck(SWIG_TYPECHECK_POINTER) Vector& {
        void* ptr = 0 ;
        if ($input && PyArray_Check($input))
            $1 = 1;
        else if (SWIG_IsOK(SWIG_ConvertPtr($input,&ptr,SWIGTYPE_p_OpenMEEG__Vector,SWIG_POINTER_EXCEPTION)))
            $1 = 1;
        else
            $1 = 0;
    }
    %typemap(in) Vector& {
        $1 = new_OpenMEEG_Vector($input);
    }
    %typemap(freearg) Vector& {
        if ($1) delete $1;
    }

    %typecheck(SWIG_TYPECHECK_POINTER) Matrix& {
        void* ptr = 0 ;
        if ($input && PyArray_Check($input))
            $1 = 1;
        else if (SWIG_IsOK(SWIG_ConvertPtr($input,&ptr,SWIGTYPE_p_OpenMEEG__Matrix,SWIG_POINTER_EXCEPTION)))
            $1 = 1;
        else
            $1 = 0;
    }
    %typemap(in) Matrix& {
        $1 = new_OpenMEEG_Matrix($input);
    }
    %typemap(freearg) Matrix& {
        if ($1) delete $1;
    }

    // C++ -> Python
    %typemap(out) unsigned& {
        $result = PyInt_FromLong(*($1));
    }
}

// /////////////////////////////////////////////////////////////////
// extensions
// /////////////////////////////////////////////////////////////////

// OpenMEEG

%ignore OpenMEEG::Filetype;

%ignore OpenMEEG::Geometry::MeshPair;  // Warning 325: Nested struct not currently supported (MeshPair ignored)
%rename(import_) import; // Warning 314: 'import' is a python keyword, renaming to '_import' (in Geometry)

%extend OpenMEEG::Geometry {

    IndexMap
    add_vertices(PyObject* pyobj) {
        return geom_add_vertices($self,pyobj);
    }
}

%extend OpenMEEG::Vertex {

    double x() { return (($self)->x()); }
    double y() { return (($self)->y()); }
    double z() { return (($self)->z()); }

    PyObject* array() const {
        npy_intp shape[1] = { 3 };
        PyObject* array = PyArray_SimpleNew(1,shape,NPY_DOUBLE);
        const double* data = static_cast<const double*>(*($self));
        double* datain = static_cast<double*>(PyArray_DATA(reinterpret_cast<PyArrayObject*>(array)));
        std::copy(&data[0],&data[3],datain);
        return array;
    }
}

%extend OpenMEEG::TriangleIndices {
    PyObject* array() const {
        npy_intp shape[1] = { 3 };
        PyObject* array = PyArray_SimpleNew(1,shape,NPY_UINT);
        const unsigned* data = &((*($self))[0]);
        unsigned* datain = static_cast<unsigned*>(PyArray_DATA(reinterpret_cast<PyArrayObject*>(array)));
        std::copy(&data[0],&data[3],datain);
        return array;
    }
}

%extend OpenMEEG::Triangle {
    double area() { return (($self)->area()); }
}

%extend OpenMEEG::Vector {
    Vector(PyObject* pyobj) { return new_OpenMEEG_Vector(pyobj); }

    PyObject* array() {
        const npy_intp ndims = 1;
        npy_intp ar_dim[] = { static_cast<npy_intp>(($self)->size()) };
        PyArrayObject* array = reinterpret_cast<PyArrayObject*>(PyArray_SimpleNewFromData(ndims,ar_dim,NPY_DOUBLE,static_cast<void*>(($self)->data())));
        return PyArray_Return(array);
    }

    // Setters

    void setvalue(const unsigned int i,const double d) { (*($self))(i) = d; }

    double value(unsigned int i) {
        if (i>=($self)->size())
            throw Error(SWIG_IndexError, "Index out of range");
        return (*($self))(i);
    }
}

%extend OpenMEEG::Matrix {

    Matrix(PyObject* pyobj) { return new_OpenMEEG_Matrix(pyobj); }

    static void Free(PyObject* capsule) {
        SharedData* shared_data = reinterpret_cast<SharedData*>(PyCapsule_GetPointer(capsule,PyCapsule_GetName(capsule)));
        delete shared_data;
    }

    PyObject* array() {

        const npy_intp ndims = 2;
        npy_intp* dims = new npy_intp[ndims];
        dims[0] = ($self)->nlin();
        dims[1] = ($self)->ncol();

        SharedData* shared_data = new SharedData(($self)->get_shared_data_ptr());
        void* data = static_cast<void*>(shared_data->get());
        PyObject* obj = PyArray_New(&PyArray_Type,ndims,dims,NPY_DOUBLE,NULL,
                                    data,0,NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_OWNDATA,NULL);

        PyArrayObject* array = reinterpret_cast<PyArrayObject*>(obj);
        if (array==nullptr)
            throw Error(SWIG_RuntimeError,"Cannot create numpy array from OpenMEEG matrix.");

        PyObject* capsule = PyCapsule_New(shared_data,"wrapped matrix",static_cast<PyCapsule_Destructor>(&OpenMEEG_Matrix_Free));
        if (PyArray_SetBaseObject(array,capsule)==-1) {
            Py_DECREF(array);
            throw Error(SWIG_RuntimeError,"Cannot create numpy array from OpenMEEG matrix.");
        }

        return PyArray_Return(array);
    }

    void setvalue(const unsigned int i,const unsigned int j,const double d) { (*($self))(i,j) = d; }

    double value(const unsigned int i,const unsigned int j) {
        if ((i>=($self)->nlin()) || (j>=($self)->ncol()))
            throw Error(SWIG_IndexError, "i or j out of range");
        return (*($self))(i,j);
    }
}

%extend OpenMEEG::SymMatrix {
    SymMatrix(PyObject* pyobj) { return new_OpenMEEG_SymMatrix(pyobj); }

    PyObject* array_flat() {
        const npy_intp ndims = 1;
        npy_intp* dims = new npy_intp[ndims];
        dims[0] = ($self)->size();

        // make a copy of the data
        double* data = new double[dims[0]];
        double* start = (*($self)).data();
        std::copy(start,start+dims[0],data);
        PyArrayObject* array = reinterpret_cast<PyArrayObject*>(PyArray_New(&PyArray_Type,ndims,dims,NPY_DOUBLE,NULL,
                                                                            static_cast<void*>(data),0,NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_OWNDATA,NULL));
        return PyArray_Return(array);
    }

    void setvalue(const unsigned int i,const unsigned int j,const double d) { (*($self))(i,j) = d; }

    double value(const unsigned int i,const unsigned int j) {
        if ((i>=($self)->nlin()) || (j>=($self)->ncol()))
            throw Error(SWIG_IndexError, "i or j out of range");
        return (*($self))(i,j);
    }
}

%ignore OpenMEEG::Mesh::name(); // ignore non const name() method

%extend OpenMEEG::Mesh {

    void add_triangles(PyObject* pyobj,const IndexMap& indmap) {
        mesh_add_triangles($self,pyobj,indmap);
    }

    Mesh(PyObject* vertices,PyObject* triangles,const std::string name="",Geometry* geom=nullptr) {
        Mesh* mesh = new Mesh(geom);
        mesh->name() = name;
        const OpenMEEG::IndexMap& indmap = geom_add_vertices(&(mesh->geometry()),vertices);
        mesh_add_triangles(mesh,triangles,indmap);
        mesh->update(true);
        return mesh;
    }

    const char* __str__() { return ($self)->name().c_str(); }
}

%extend OpenMEEG::Geometry {

    Geometry(PyObject* pylist) {

        if (pylist==nullptr || !PyList_Check(pylist))
            throw Error(SWIG_TypeError, "Argument to Geometry constructor must be a list");

        //  Add vertices of all meshes.

        const unsigned N = PyList_Size(pylist);
        if (N==0)
            throw Error(SWIG_ValueError, "Argument to Geometry constructor must be a non-empty list");
        OpenMEEG::Geometry* geometry = new OpenMEEG::Geometry(N);

        std::vector<OpenMEEG::IndexMap> indmap(N);
        for (unsigned i=0; i<N; ++i) {
            PyObject* item = PyList_GetItem(pylist,i);
            if (item==nullptr || !PyList_Check(item) || PyList_Size(item)!=3)
                throw Error(SWIG_TypeError, "Geometry constructor argument must be a list of lists, each of length 3");
            PyObject* vertices = PyList_GetItem(item,1);
            indmap[i] = geom_add_vertices(geometry,vertices);
        }

        //  Create meshes and add triangles.

        for (unsigned i=0; i<N; ++i) {
            PyObject* item = PyList_GetItem(pylist,i);
            PyObject* name = PyList_GetItem(item,0);
            if (name==nullptr || !PyUnicode_Check(name))
                throw Error(SWIG_TypeError, "Geometry constructor list of lists must each have first entry a non-empty string.");
            Mesh& mesh = geometry->add_mesh(PyUnicode_AsUTF8AndSize(name, nullptr));  // Since 3.10
            PyObject* triangles = PyList_GetItem(item,2);
            mesh_add_triangles(&mesh,triangles,indmap[i]);
            mesh.update(true);
            #ifdef DEBUG
            std::ostringstream oss;
            oss << "SWIG Geometry update of meshes(" << i << ")=\"" << mesh.name() << "\"";
            for (const auto& mesh : geometry->meshes())
                mesh.check_consistency(oss.str());
            #endif
        }

        return geometry;
    }
}

// Input

%include <logger.h>
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
%include <GeometryIO.h>
%include <sensors.h>
%include <mesh.h>
%include <integrator.h>
%include <interface.h>
%include <domain.h>
%include <assemble.h>
%include <gain.h>
%include <forward.h>

#ifdef _MSC_VER
#pragma warning( pop )
#endif
