#ifndef H_MESH
#define H_MESH

#include <istream>
#include <fstream>
#include <set>
#include "vect3.h"
#include "triangle.h"

#ifdef USE_VTK
#include <vtkPolyDataReader.h>
#endif

typedef std::set<int> intSet;

/** \brief  Mesh

    Mesh Class

© Copyright 2007-2007 Odyssee INRIA . All Rights Reserved.

    \author $LastChangedBy$
    \date $LastChangedDate$
    \version $Rev$  \sa
**/

class Mesh
{

private:
    int npts; //!< Number of points in the mesh.
    int ntrgs; //!< Number of triangles in the mesh.
    Vect3* pts; //!< Points of the mesh
    Triangle* trgs; //!< Triangles of the mesh
    intSet* links; //!< links[i] are the triangles that contain point i : each point knows each triangle it is a part of
    Vect3* normals; //!< Normals at each point

public:

    enum Filetype { VTK, TRI, BND };

    Mesh();
    Mesh(int, int); // npts, ntrgs
    Mesh(int, Vect3 *, int ,Triangle* );
    Mesh(const Mesh& M);
    //Mesh& operator=( Mesh& M);
    Mesh& operator=(const Mesh& M);
    ~Mesh(){kill();}
    inline int nbPts() const { return npts; }
    inline int nbTrgs() const { return ntrgs; }
    inline const Vect3& getPt(int i) const { return pts[i]; }
    inline const Triangle& getTrg(int i) const { return trgs[i]; }
    inline const intSet& getTrianglesForPoint(int i) const { return links[i]; }
    inline intSet* getTrianglesForPoints() const { return links; }

    inline Vect3 normal(int i) const { return normals[i]; }
    inline Vect3& normal(int i) { return normals[i]; }

    inline Vect3& operator[](int i) {
        assert(i<npts);
        return pts[i];
    }

    inline void setPt(int i,const Vect3& v) { pts[i]=v;}
    inline void setTrg(int i,const Triangle& t) { trgs[i]=t;}

    void make_links();
    void compute_omega();

    /**
       * Get center of triangle
       * \param i index of the triangle
       */
    Vect3 center(int) const;

    /**
       * Read mesh from file
       * \param filename can be .vtk, .tri (ascii), .bnd
       * \param verbose true or false
       */
    void load(const char* filename, bool verbose=true);

#ifdef USE_VTK
    void load_vtk(std::istream &is);
    void load_vtk(const char*);
#else
    template <typename T>
    void load_vtk(T) {
        std::cerr << "You have to compile OpenMEEG with VTK to read VTK files" << std::endl;
        exit(1);
    }
#endif

    void load_tri(std::istream &);
    void load_tri(const char*);

    void load_bnd(std::istream &);
    void load_bnd(const char*);

    void load_Mesh(std::istream &is);
    void load_Mesh(const char*);

    void save(const char* filename);
    void save_vtk(const char*);
    void save_bnd(const char*);
    void save_tri(const char*);
    void save_mesh(const char*);

    /** \brief elem

            Find the triangles that use the point i
            and store the indices of these points in T

        \param i point ID
        \param T indices of the triangles that use the point i
        \return the number of triangles that use point i
        \sa
    **/
    inline int elem(int i, int* T ) const{ // gives each triangle a point belongs to
        int c=0;
        for(int k=0; k<ntrgs; k++){
            if( (i==trgs[k].s1()) || (i==trgs[k].s2()) || (i==trgs[k].s3())){
                T[c] = k;
                c = c+1;
            }
        }
        return c;
    }

    void getFileFormat(const char* filename);

    /** \brief append a mesh

            append a mesh to the current Mesh.
            Used to concat two meshes

        \param m Mesh to append
        \return void
        \sa
    **/
    void append(const Mesh* m);

    inline friend void operator>>(std::istream &ifs,Mesh &m){
        Filetype format = m.streamFormat;
        switch (format){
            case TRI:       m.load_tri(ifs); break;
#ifdef USE_VTK
            case VTK:       m.load_vtk(ifs); break;
#endif
            case BND:       m.load_bnd(ifs); break;
            default: std::cout << "Unknown file format" << std::endl;
        }
    }

private:
    Filetype streamFormat;
    void kill();
    void copy(const Mesh& M);
#ifdef USE_VTK
    void getDataFromVTKReader(vtkPolyDataReader* vtkMesh);
#endif
};

#endif
