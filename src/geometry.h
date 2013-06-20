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

#ifndef OPENMEEG_GEOMETRY_H
#define OPENMEEG_GEOMETRY_H

#include <assert.h>
#include <set>
#include <vector>
#include "vertex.h"
#include "triangle.h"
#include "mesh.h"
#include "interface.h"
#include "domain.h"

// for IO:s
#include <iostream>
#include <fstream>
#include <vector>
#include <new>
#include <string>
#ifdef USE_VTK
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkDataReader.h>
#include <vtkCellArray.h>
#include <vtkCharArray.h>
#include <vtkProperty.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkCleanPolyData.h>
#endif

namespace OpenMEEG {

    /** \brief  Geometry
        Geometry Class
    **/

    enum Filetype { VTP, VTK, TRI, BND, MESH, OFF, GIFTI };


    class OPENMEEG_EXPORT Geometry {

        // Geometry types
        typedef std::vector<Vect3>       Normals;
        // typedef std::vector<std::string> Names;

    public:

        friend class Reader;

        // Default iterator of a Geometry is an Iterator on the meshes
        typedef Meshes::iterator       iterator;
        typedef Meshes::const_iterator const_iterator;

        // Iterators
              iterator             begin()                 { return meshes().begin();    }
        const_iterator             begin()           const { return meshes().begin();    }
              iterator             end()                   { return meshes().end();      }
        const_iterator             end()             const { return meshes().end();      }
        Vertices::iterator         vertex_begin()          { return vertices().begin();  }
        Vertices::const_iterator   vertex_begin()    const { return vertices().begin();  }
        Vertices::iterator         vertex_end()            { return vertices().end();    }
        Vertices::const_iterator   vertex_end()      const { return vertices().end();    }
        Domains::iterator          domain_begin()          { return domains().begin();   }
        Domains::const_iterator    domain_begin()    const { return domains().begin();   }
        Domains::iterator          domain_end()            { return domains().end();     }
        Domains::const_iterator    domain_end()      const { return domains().end();     }
        Interfaces::iterator       interface_begin()       { return interfaces().begin();}
        Interfaces::const_iterator interface_begin() const { return interfaces().begin();}
        Interfaces::iterator       interface_end()         { return interfaces().end();  }
        Interfaces::const_iterator interface_end()   const { return interfaces().end();  }

        Geometry():  _has_cond(false)   { }
        ~Geometry()                     { destroy(); }
        void          read(const char* geomFileName, const char* condFileName = NULL);

              bool&         has_cond()                    { return _has_cond; }
        const bool&         has_cond()              const { return _has_cond; }
              Vertices&     vertices()                    { return _vertices; }
        const Vertices&     vertices()              const { return _vertices; }
              Normals&      normals()                     { return _normals; }
        const Normals&      normals()               const { return _normals; }
              Meshes&       meshes()                      { return _meshes; }
        const Meshes&       meshes()                const { return _meshes; }
              Interfaces&   interfaces()                  { return _interfaces; }
        const Interfaces&   interfaces()            const { return _interfaces; }
              Domains&      domains()                     { return _domains; }
        const Domains&      domains()               const { return _domains; }
              size_t&       size()                        { return _size; }
        const size_t        size()                  const { return _size; }
        const size_t        nb_domains()            const { return _domains.size(); }
        
        inline double sigma      (const Domain d)                 const { return (d.sigma()); }
        inline double sigma      (const Mesh& m1, const Mesh& m2) const { return funct_on_domains(m1, m2, '+');} // return the (sum) conductivity(ies) of the shared domain(s).
        inline double sigma_diff (const Mesh& m1, const Mesh& m2) const { return funct_on_domains(m1, m2, '-');} // return the (difference) of conductivity(ies) of the shared domain(s).
        inline double sigma_inv  (const Mesh& m1, const Mesh& m2) const { return funct_on_domains(m1, m2, '/');} // return the (sum) inverse of conductivity(ies) of the shared domain(s).
        inline double indicatrice(const Mesh& m1, const Mesh& m2) const { return funct_on_domains(m1, m2, '1');} // return the (sum) indicatrice function of the shared domain(s).
        inline double sigma      (const std::string&) const;
        // bool selfCheck() const;
        // bool check(const Mesh& m) const;
        const Domain  get_domain(const Vect3&) const;
        double oriented(const Mesh&, const Mesh&) const;

        // for IO:s
        void load_tri(std::istream &, Mesh &);
        void load_tri(const char*, Mesh &);
        void load_bnd(std::istream &, Mesh &);
        void load_bnd(const char*, Mesh &);
        void load_off(std::istream &, Mesh &);
        void load_off(const char*, Mesh &);
        void load_mesh(std::istream &, Mesh &);
        void load_mesh(const char*, Mesh &);
        #ifdef USE_VTK
        void load_vtp(std::istream &);
        void load_vtp(const char*);
        void load_vtk(std::istream &, Mesh &);
        void load_vtk(const char*, Mesh &);
        void get_data_from_vtk_reader(vtkPolyDataReader* vtkMesh, Mesh &);
        #else
        template <typename T>
        void load_vtp(T) {
            std::cerr << "You have to compile OpenMEEG with VTK to read VTK/VTP files" << std::endl;
            exit(1);
        }
        template <typename T>
        void load_vtk(T, Mesh &) {
            std::cerr << "You have to compile OpenMEEG with VTK to read VTK/VTP files" << std::endl;
            exit(1);
        }
        #endif
        #ifdef USE_GIFTI
        void load_gifti(std::istream &, Mesh &);
        void load_gifti(const char*, Mesh &);
        void save_gifti(const char*, Mesh &);
        gifti_image* to_gifti_image();
        void from_gifti_image(gifti_image* gim);
        #else
        template <typename T>
        void load_gifti(T, Mesh &) {
            std::cerr << "You have to compile OpenMEEG with GIFTI to read GIFTI files" << std::endl;
            exit(1);
        }
        template <typename T>
        void save_gifti(T, Mesh &) {
            std::cerr << "You have to compile OpenMEEG with GIFTI to save GIFTI files" << std::endl;
            exit(1);
        }
        #endif
        void save_vtk(const char*)  const;
        void save_bnd(const char*)  const;
        void save_tri(const char*)  const;
        void save_off(const char*)  const;
        void save_mesh(const char*) const;
        void mesh_load(const char* filename, Mesh &m) {
            std::string extension = getNameExtension(filename);
            std::transform(extension.begin(), extension.end(), extension.begin(), (int(*)(int))std::tolower);
            if (extension == std::string("vtk")) {
                load_vtk(filename, m);
            } else if (extension == std::string("vtp")) {
                load_vtp(filename);
            } else if (extension == std::string("tri")) {
                load_tri(filename, m);
            } else if (extension == std::string("bnd")) {
                load_bnd(filename, m);
            } else if (extension == std::string("mesh")) {
                load_mesh(filename, m);
            } else if (extension == std::string("off")) {
                load_off(filename, m);
            } else if (extension == std::string("gii")) {
                load_gifti(filename, m);
            } else {
                std::cerr << "IO: load: Unknown mesh file format for " << filename << std::endl;
                exit(1);
            }
        }

    private:

        // Members
        Vertices   _vertices;
        Normals    _normals;
        Meshes     _meshes;
        Interfaces _interfaces;
        Domains    _domains;
        size_t     _size;   // total number = nb of vertices + nb of triangles
        bool       _has_cond;

        void read_geom(const std::string&);
        void read_cond(const char*);
        void geom_generate_indices();

        bool is_relative_path(const std::string& name);
        #if WIN32
        static const char PathSeparator[] = "/\\";
        #else
        static const char PathSeparator   = '/';
        #endif

        inline Domains common_domains(const Mesh& m1, const Mesh& m2) const {
            // std::set<Domain&> sdom1;
            // std::set<Domain&> sdom2;
            // Domains doms;
            // Domains::iterator dit;
            // for (dit = domains.begin(); dit != domains.end(); dit++) {
                // if (dit->meshOrient(m1) != 0) {
                    // sdom1.insert(std::abs(dit->second));
                // }
                // if (dit->meshOrient(m2) != 0) {
                    // sdom2.insert(std::abs(dit->second));
                // }
            // }
            // dit = std::set_intersection(sdom1.begin(), sdom1.end(), sdom2.begin(), sdom2.end(), doms );
            // return doms;
        }

        inline double funct_on_domains(const Mesh& m1, const Mesh& m2, char f) const {       // return the (sum) conductivity(ies) of the shared domain(s).
            Domains doms = common_domains(m1, m2);
            double ans=0.;
            for (Domains::iterator dit = doms.begin(); dit != doms.end(); dit++) {
                if (f=='+') {
                    ans += dit->sigma();
                } else if (f == '-') {
                    ans -= -1.*dit->sigma();
                } else if (f == '/') {
                    ans += 1./dit->sigma();
                } else {
                    ans += 1.;
                }
            }
            return ans;
        }

        void destroy() {
            // if (n!=0) {
                // n=0;
                // delete []M;
                // if(has_cond)
                // {
                    // delete []sigin;
                    // delete []sigout;
                // }
            // }
        }
    };

}

#endif  //! OPENMEEG_GEOMETRY_H
