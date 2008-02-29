/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
*/

#ifndef H_MESH
#define H_MESH

#include <istream>
#include <fstream>
#include <set>
#include "vect3.h"
#include "triangle.h"
#include "sparse_matrice.h"
#include "symmatrice.h"

#ifdef USE_VTK
#include <vtkPolyDataReader.h>
#endif

typedef std::set<int> intSet;
int tri_tri_overlap_test_3d(float p1[3], float q1[3], float r1[3],float p2[3], float q2[3], float r2[3]);

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
    int ncomponents; //!< Number of connexe components
    int ninversions; //!< Number of triangles with inverted orientations

    void load_tri(std::istream &, bool checkOrientations = true);
    void load_tri(const char*, bool checkOrientations = true);

    void load_bnd(std::istream &, bool checkOrientations = true);
    void load_bnd(const char*, bool checkOrientations = true);

    void load_mesh(std::istream &is, bool checkOrientations = true);
    void load_mesh(const char*, bool checkOrientations = true);
    #ifdef USE_VTK
        void load_vtk(std::istream &is, bool checkOrientations = true);
        void load_vtk(const char*, bool checkOrientations = true);
    #else
        template <typename T>
        void load_vtk(T, bool checkOrientations = true) {
            std::cerr << "You have to compile OpenMEEG with VTK to read VTK files" << std::endl;
            exit(1);
        }
    #endif

    void save_vtk(const char*);
    void save_bnd(const char*);
    void save_tri(const char*);
    void save_mesh(const char*);

    void updateTriangleOrientations(bool checkOrientations = true);
    void update_triangles();
    void recompute_normals();

public:

    enum Filetype { VTK, TRI, BND, MESH };

    Mesh();
    Mesh(int, int); // npts, ntrgs
    Mesh(const Mesh& M);
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
    // void compute_omega();

    /**
       * Get center of triangle
       * \param i index of the triangle
       */
    // Vect3 center(int) const;

    /**
       * Read mesh from file
       * \param filename can be .vtk, .tri (ascii), .bnd
       * \param checkOrientations true or false
       * \param verbose true or false
       */
    void load(const char* filename, bool checkOrientations = true, bool verbose=true);

    /**
       * Save mesh to file
       * \param filename can be .vtk, .tri (ascii), .bnd or .mesh
       */
    void save(const char* filename);

    /** \brief Get file format based on file extension

            A list of supported file formats is in variable "Filetype"

        \param filename
        \return void
        \sa
    **/
    void getFileFormat(const char* filename);

    /** \brief append a mesh

            append a mesh to the current Mesh.
            Used to concat two meshes

        \param m Mesh to append
        \return void
        \sa
    **/
    void append(const Mesh* m);
    
    /** \brief Print info

            Print to std::cout some info about the mesh

        \param m Mesh to append
        \return void
        \sa
    **/
    void info() const;
    
    /** \brief Smooth Mesh

            Smooth Mesh

        \param smoothing_intensity
        \param niter
        \return void
        \sa
    **/
    void smooth(double smoothing_intensity, size_t niter);
    
    int getNeighTrg(int a, int b, int c) const;
    
    /**
     * Compute the surfacic gradient
    **/
    sparse_matrice gradient() const;

    /**
     * Return the area of each triangle
    **/
    vecteur areas() const;

    inline friend void operator>>(std::istream &ifs,Mesh &m){
        Filetype format = m.streamFormat;
        switch (format){
            case TRI:       m.load_tri(ifs); break;
#ifdef USE_VTK
            case VTK:       m.load_vtk(ifs); break;
#endif
            case BND:       m.load_bnd(ifs); break;
            case MESH:      m.load_mesh(ifs); break;
            default: std::cout << "Unknown file format" << std::endl;
        }
    }
    
    /**
     * Check intersection between two triangles
    **/
    static inline bool triangle_intersection( const Mesh& m1, int nT1, const Mesh& m2, int nT2 ) {
        const Triangle& T1 = m1.getTrg(nT1);
        const Vect3& p1 = m1.getPt(T1.s1());
        const Vect3& q1 = m1.getPt(T1.s2());
        const Vect3& r1 = m1.getPt(T1.s3());
        const Triangle& T2 = m2.getTrg(nT2);
        const Vect3& p2 = m2.getPt(T2.s1());
        const Vect3& q2 = m2.getPt(T2.s2());
        const Vect3& r2 = m2.getPt(T2.s3());

        float pp1[3] = {p1.x(),p1.y(),p1.z()};
        float qq1[3] = {q1.x(),q1.y(),q1.z()};
        float rr1[3] = {r1.x(),r1.y(),r1.z()};
        float pp2[3] = {p2.x(),p2.y(),p2.z()};
        float qq2[3] = {q2.x(),q2.y(),q2.z()};
        float rr2[3] = {r2.x(),r2.y(),r2.z()};
        return tri_tri_overlap_test_3d(pp1,qq1,rr1,pp2,qq2,rr2);
    }

    bool selfIntersection() const;
    bool intersection(const Mesh& m) const;

private:
    Filetype streamFormat;
    void kill();
    void copy(const Mesh& M);
#ifdef USE_VTK
    void getDataFromVTKReader(vtkPolyDataReader* vtkMesh);
#endif
};

/**
 * P1Vector : aux function to compute the surfacic gradient
**/
inline Vect3 P1Vector( const Vect3 &p0, const Vect3 &p1, const Vect3 &p2, const int idx )
{
    assert(idx>-1 && idx<3);
    int i = idx+1;
    Vect3 pts[5] = {p2,p0,p1,p2,p0};
    Vect3 ret(0,0,0);
    Vect3 pim1pi = pts[i]-pts[i-1];
    Vect3 pim1pip1 = pts[i+1]-pts[i-1];
    Vect3 pim1H = ( (1.0/pim1pip1.norme2()) * ( pim1pi*pim1pip1 ) ) * pim1pip1;
    Vect3 piH = pim1H-pim1pi;
    ret = -1.0/piH.norme2()*piH;

    return ret;
}

#endif
