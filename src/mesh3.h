/*
   Project Name : OpenMEEG

   © INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
   GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
   Emmanuel OLIVI
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

#ifdef USE_GIFTI
extern "C" {
#include <gifti_io.h>
}
#endif

namespace OpenMEEG
{

    typedef std::set<int> intSet;
    OPENMEEG_EXPORT bool tri_tri_overlap_test_3d(double p1[3], double q1[3], double r1[3],
            double p2[3], double q2[3], double r2[3]);


    class OPENMEEG_EXPORT Mesh
    {

        public:

            void update_triangles();
            void recompute_normals();
            void flip_faces();


            enum Filetype { VTK, TRI, BND, MESH, OFF, GIFTI };

            Mesh();
            Mesh(const int, const int); // npts, ntrgs
            Mesh(const Mesh& M);
            Mesh& operator=(const Mesh& M);

            ~Mesh() { destroy(); }

            int nbPts()  const { return npts;       }
            int nbTrgs() const { return ntrgs;      }
            int size()   const { return npts+ntrgs; }

            const Vect3&    point(const int i)  const    { return pts[i];  }
                  Vect3&    point(const int i)           { return pts[i];  }
            const Triangle& triangle(const int i) const  { return trgs[i]; }
            Triangle& triangle(const int i)              { return trgs[i]; }

            const intSet& get_triangles_for_point(const int i) const { return links[i]; }
            intSet* get_triangles_for_points()           const { return links;    }

            Vect3  normal(const int i) const { return normals[i]; }
            Vect3& normal(const int i)       { return normals[i]; }

            Vect3& operator[](int i) {
                assert(i<npts);
                return pts[i];
            }

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
            void save(const char* filename) const;

            /** \brief Get file format based on file extension

              A list of supported file formats is in variable "Filetype"

              \param filename
              \return void
              \sa
             **/
            void get_file_format(const char* filename);

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

            int neighbor_triangle(int a, int b, int c) const;

            /**
             * Compute the surfacic gradient
             **/
            SparseMatrix gradient() const;

            /**
             * Compute the total solid angle of a surface for a point v and tell whether v is inside the mesh or not.
             **/
            // bool contains_point(const Vect3&) const;

            /**
             * Return the area of each triangle
             **/
            Vector areas() const;

            /**
             * Check intersection between two triangles
             **/
            static inline bool triangle_intersection( const Mesh& m1, int nT1, const Mesh& m2, int nT2 ) {
                const Triangle& T1 = m1.triangle(nT1);
                const Vect3& p1 = m1.point(T1.s1());
                const Vect3& q1 = m1.point(T1.s2());
                const Vect3& r1 = m1.point(T1.s3());
                const Triangle& T2 = m2.triangle(nT2);
                const Vect3& p2 = m2.point(T2.s1());
                const Vect3& q2 = m2.point(T2.s2());
                const Vect3& r2 = m2.point(T2.s3());

                double pp1[3] = {p1.x(), p1.y(), p1.z()};
                double qq1[3] = {q1.x(), q1.y(), q1.z()};
                double rr1[3] = {r1.x(), r1.y(), r1.z()};
                double pp2[3] = {p2.x(), p2.y(), p2.z()};
                double qq2[3] = {q2.x(), q2.y(), q2.z()};
                double rr2[3] = {r2.x(), r2.y(), r2.z()};
                return tri_tri_overlap_test_3d(pp1, qq1, rr1, pp2, qq2, rr2);
            }

            bool has_self_intersection() const;
            bool intersection(const Mesh& m) const;
            bool has_correct_orientation() const;

        private:
            int npts; //!< Number of points in the mesh.
            int ntrgs; //!< Number of triangles in the mesh.
            Vect3* pts; //!< Points of the mesh
            Triangle* trgs; //!< Triangles of the mesh
            intSet* links; //!< links[i] are the triangles that contain point i : each point knows each triangle it is a part of
            Vect3* normals; //!< Normals at each point

    };

    /**
     * P1Vector : aux function to compute the surfacic gradient // TODO a quoi ça sert ?
     **/
    inline Vect3 P1Vector( const Vect3 &p0, const Vect3 &p1, const Vect3 &p2, const int idx )
    {
        assert(idx>-1 && idx<3);
        int i = idx+1;
        Vect3 pts[5] = {p2, p0, p1, p2, p0};
        Vect3 ret(0, 0, 0);
        Vect3 pim1pi = pts[i]-pts[i-1];
        Vect3 pim1pip1 = pts[i+1]-pts[i-1];
        Vect3 pim1H = ( (1.0/pim1pip1.norm2()) * ( pim1pi*pim1pip1 ) ) * pim1pip1;
        Vect3 piH = pim1H-pim1pi;
        ret = -1.0/piH.norm2()*piH;

        return ret;
    }

} // end namespace OpenMEEG

#endif  //! OPENMEEG_MESH_H
