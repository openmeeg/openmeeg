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

// for IO:s
#include <iostream>
#include <fstream>

#include <vector>
#include <set>
#include <string>
#include "triangle.h"
#include <IOUtils.H>
#include "om_utils.h"

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

    /** \brief  Mesh

        Mesh class

    **/

    class OPENMEEG_EXPORT Mesh: public Triangles {

    public:
        
        typedef Triangles                                   base;
        typedef std::set<Vertex *, std::less<Vertex *> >    SetPVertex;
        typedef std::vector<Vertex *>                       VectPVertex;
        typedef VectPVertex::iterator                       vertex_iterator;
        typedef VectPVertex::const_iterator                 const_vertex_iterator;

        Mesh(): base(), _name(""), _all_vertices(NULL), _outermost(false) {};

        Mesh(Vertices &all_vertices, Normals &all_normals, const std::string name = ""): _all_vertices(&all_vertices), _all_normals(&all_normals), _name(name), _outermost(false) { }
        Mesh(Vertices &all_vertices, const std::string name = ""): _all_vertices(&all_vertices), _name(name), _outermost(false) { }

        Mesh(std::string name): _name(name), _all_vertices(NULL), _outermost(false) { }

        void recompute_normals();

        // Iterators on vertices
              vertex_iterator      vertex_begin()                { return vertices().begin(); }
        const_vertex_iterator      vertex_begin()  const         { return vertices().begin(); }
              vertex_iterator      vertex_end()                  { return vertices().end(); }
        const_vertex_iterator      vertex_end()    const         { return vertices().end(); }

        std::vector<SetTriangle>   links()         const         { return _links; }
        std::vector<SetTriangle> & links()                       { return _links; }

        std::string                name()          const         { return _name; }
        std::string &              name()                        { return _name; }

        void                       add_vertex(const Vertex &v)   { _all_vertices->push_back(v); _vertices.push_back(&*_all_vertices->end());}
        void                       add_normal(const Normal &n)   { _all_normals->push_back(n);}

        VectPVertex                vertices()      const         { return _vertices; }
        VectPVertex &              vertices()                    { return _vertices; }

        const size_t               nb_vertices()   const         { return _vertices.size(); }
        const size_t               nb_triangles()  const         { return this->size(); }

          Normals                all_normals()  const         { return *_all_normals; }
          Normals                all_normals()                { return *_all_normals; }
          Vertices                 all_vertices()  const         { return *_all_vertices; }
          Vertices                 all_vertices()                { return *_all_vertices; }

        // mesh state
        void info() const ;
        bool triangle_intersection(const Triangle&, const Triangle&) const;
        bool has_self_intersection() const;
        bool intersection(const Mesh&) const;
        bool has_correct_orientation() const;
        void update();

        const SetTriangle& get_triangles_for_point(const Vertex& V) const {
            size_t i = 0;
            for (const_vertex_iterator vit = vertices().begin(); vit != vertices().end(); vit++, i++) {
                if (*vit == &V) {
                    return links()[i];
                }
            }
        }

        friend std::istream& operator>>(std::istream &is, Mesh &m);

        //  Returns True if it is an outermost mesh.
              bool&        outermost()       { return _outermost; }
        const bool&        outermost() const { return _outermost; }

        // for IO:s
        void load_mesh(const char* filename, const bool &verbose = true);
        void load_tri_file(std::istream &);
        void load_tri_file(const char*);
        void load_bnd_file(std::istream &);
        void load_bnd_file(const char*);
        void load_off_file(std::istream &);
        void load_off_file(const char*);
        void load_mesh_file(std::istream &);
        void load_mesh_file(const char*);
        #ifdef USE_VTK
        void load_vtp_file(std::istream &);
        void load_vtp_file(const char*);
        void load_vtk_file(std::istream &);
        void load_vtk_file(const char*);
        void get_data_from_vtk_reader(vtkPolyDataReader* vtkMesh);
        #else
        template <typename T>
        void load_vtp_file(T) {
            std::cerr << "You have to compile OpenMEEG with VTK to read VTK/VTP files" << std::endl;
            exit(1);
        }
        template <typename T>
        void load_vtk_file(T) {
            std::cerr << "You have to compile OpenMEEG with VTK to read VTK/VTP files" << std::endl;
            exit(1);
        }
        #endif
        #ifdef USE_GIFTI
        void load_gifti_file(std::istream &);
        void load_gifti_file(const char*);
        void save_gifti_file(const char*);
        gifti_image* to_gifti_image();
        void from_gifti_image(gifti_image* gim);
        #else
        template <typename T>
        void load_gifti_file(T) {
            std::cerr << "You have to compile OpenMEEG with GIFTI to read GIFTI files" << std::endl;
            exit(1);
        }
        template <typename T>
        void save_gifti_file(T) {
            std::cerr << "You have to compile OpenMEEG with GIFTI to save GIFTI files" << std::endl;
            exit(1);
        }
        #endif
        void save_vtk_file(const char*)  const;
        void save_bnd_file(const char*)  const;
        void save_tri_file(const char*)  const;
        void save_off_file(const char*)  const;
        void save_mesh_file(const char*) const;
        void save_mesh(const char*) const ;

        void destroy() {}; // TODO
    private:
        
        std::string                 _name;
        std::vector<SetTriangle>    _links;
        Vertices *                  _all_vertices;
        VectPVertex                 _vertices;
        Normals *                _all_normals;
        bool                        _outermost;

    };

    typedef std::vector<Mesh>        Meshes;
}

#endif  //  ! OPENMEEG_MESH_H
