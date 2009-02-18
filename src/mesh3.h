/* FILE: $Id$ */

/*
Project Name : OpenMEEG

version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Maureen.Clerc.AT.sophia.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.sophia.inria.fr)

The OpenMEEG software is a C++ package for solving the forward/inverse
problems of electroencephalography and magnetoencephalography.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's authors,  the holders of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

#ifndef OPENMEEG_MESH_H
#define OPENMEEG_MESH_H

#include <istream>
#include <fstream>
#include <set>
#include "vect3.h"
#include "triangle.h"
#include "vector.h"
#include "matrix.h"
#include "symmatrix.h"
#include "sparse_matrix.h"

#ifdef USE_VTK
#include <vtkPolyDataReader.h>
#endif

namespace OpenMEEG {

    typedef std::set<int> intSet;
    bool tri_tri_overlap_test_3d(double p1[3], double q1[3], double r1[3],double p2[3], double q2[3], double r2[3]);


    class OPENMEEG_EXPORT Mesh
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

        void load_tri(std::istream &);
        void load_tri(const char*);

        void load_bnd(std::istream &);
        void load_bnd(const char*);

        void load_off(std::istream &);
        void load_off(const char*);

        void load_mesh(std::istream &is);
        void load_mesh(const char*);

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

        void save_vtk(const char*);
        void save_bnd(const char*);
        void save_tri(const char*);
        void save_off(const char*);
        void save_mesh(const char*);

        void update_triangles();
        void recompute_normals();

    public:

        enum Filetype { VTK, TRI, BND, MESH, OFF };

        Mesh();
        Mesh(int, int); // npts, ntrgs
        Mesh(const Mesh& M);
        Mesh& operator=(const Mesh& M);
        ~Mesh(){kill();}
        inline int nbPts() const { return npts; }
        inline int nbTrgs() const { return ntrgs; }
        inline const Vect3& getPt(int i) const { return pts[i]; }
        inline const Triangle& getTrg(int i) const { return trgs[i]; }
        inline Triangle& getTrg(int i) { return trgs[i]; }
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

        /**
           * Get center of triangle
           * \param i index of the triangle
           */
        // Vect3 center(int) const;

        /**
           * Read mesh from file
           * \param filename can be .vtk, .tri (ascii), .bnd
           * \param verbose true or false
           */
        void load(const char* filename, bool verbose=true);

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
        SparseMatrix gradient() const;

        /**
         * Return the area of each triangle
        **/
        Vector areas() const;

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

            double pp1[3] = {p1.x(),p1.y(),p1.z()};
            double qq1[3] = {q1.x(),q1.y(),q1.z()};
            double rr1[3] = {r1.x(),r1.y(),r1.z()};
            double pp2[3] = {p2.x(),p2.y(),p2.z()};
            double qq2[3] = {q2.x(),q2.y(),q2.z()};
            double rr2[3] = {r2.x(),r2.y(),r2.z()};
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
        Vect3 pim1H = ( (1.0/pim1pip1.norm2()) * ( pim1pi*pim1pip1 ) ) * pim1pip1;
        Vect3 piH = pim1H-pim1pi;
        ret = -1.0/piH.norm2()*piH;

        return ret;
    }
}

#endif  //! OPENMEEG_MESH_H
