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

    enum Filetype { VTK, TRI, BND, MESH, OFF, GIFTI };

    class OPENMEEG_EXPORT Mesh: public Triangles {

    public:
        
        typedef std::vector<Vertex *>                       VectPVertex;
        typedef VectPVertex::const_iterator                 const_vertex_iterator;
        typedef VectPVertex::const_reverse_iterator         const_vertex_reverse_iterator;

        // Constructors
        Mesh(): Triangles(), name_(""), all_vertices_(NULL), outermost_(false) {};

        Mesh(Vertices &all_vertices, const std::string name = ""): all_vertices_(&all_vertices), name_(name), outermost_(false) { }

        Mesh(std::string name): name_(name), all_vertices_(NULL), outermost_(false) { }

        // Destructor
        ~Mesh() { }

        // Iterators on vertices
        const_vertex_reverse_iterator vertex_rbegin() const         { return vertices_.rbegin(); }
        const_vertex_reverse_iterator vertex_rend()   const         { return vertices_.rend(); }
        const_vertex_iterator         vertex_begin()  const         { return vertices_.begin(); }
        const_vertex_iterator         vertex_end()    const         { return vertices_.end(); }

        std::vector<SetTriangle>&     links()                       { return links_; }
        std::vector<SetTriangle>      links()         const         { return links_; }

        std::string                   name()          const         { return name_; }
        std::string &                 name()                        { return name_; }

        void                          add_vertex(const Vertex &v)   { all_vertices_->push_back(v); vertices_.push_back(&(*all_vertices_->rbegin())); }

        VectPVertex                   vertices()      const         { return vertices_; }
        VectPVertex &                 vertices()                    { return vertices_; }

        const size_t                  nb_vertices()   const         { return vertices_.size(); }
        const size_t                  nb_triangles()  const         { return this->size(); }

              Vertices                all_vertices()  const         { return *all_vertices_; }

        // Mesh state
        void info() const ;
        bool has_self_intersection() const;
        bool intersection(const Mesh&) const;
        bool has_correct_orientation() const;
        bool triangle_intersection(const Triangle&, const Triangle&) const;
        void update();

        const SetTriangle& get_triangles_for_point(const Vertex& V) const;

        friend std::istream& operator>>(std::istream &is, Mesh &m);

        //  Returns True if it is an outermost mesh.
              bool&        outermost()       { return outermost_; }
        const bool&        outermost() const { return outermost_; }

        // for IO:s
        size_t load_mesh(const char* filename, const bool &verbose = true, const bool &read_all = true);
        size_t load_tri_file(std::istream &, const bool &read_all = true);
        size_t load_tri_file(const char*, const bool &read_all = true);
        size_t load_bnd_file(std::istream &, const bool &read_all = true);
        size_t load_bnd_file(const char*, const bool &read_all = true);
        size_t load_off_file(std::istream &, const bool &read_all = true);
        size_t load_off_file(const char*, const bool &read_all = true);
        size_t load_mesh_file(std::istream &, const bool &read_all = true);
        size_t load_mesh_file(const char*, const bool &read_all = true);
        #ifdef USE_VTK
        size_t load_vtp_file(std::istream &, const bool &read_all = true);
        size_t load_vtp_file(const char*, const bool &read_all = true);
        size_t load_vtk_file(std::istream &, const bool &read_all = true);
        size_t load_vtk_file(const char*, const bool &read_all = true);
        size_t get_data_from_vtk_reader(vtkPolyDataReader* vtkMesh);
        #else
        template <typename T>
        size_t load_vtp_file(T, const bool &read_all = true) {
            std::cerr << "You have to compile OpenMEEG with VTK to read VTK/VTP files" << std::endl;
            exit(1);
        }
        template <typename T>
        size_t load_vtk_file(T, const bool &read_all = true) {
            std::cerr << "You have to compile OpenMEEG with VTK to read VTK/VTP files" << std::endl;
            exit(1);
        }
        #endif
        #ifdef USE_GIFTI
        size_t load_gifti_file(std::istream &, const bool &read_all = true);
        size_t load_gifti_file(const char*, const bool &read_all = true);
        void save_gifti_file(const char*);
        gifti_image* to_gifti_image();
        void from_gifti_image(gifti_image* gim);
        #else
        template <typename T>
        size_t load_gifti_file(T, const bool &read_all = true) {
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

    private:
        
        std::string                 name_;
        std::vector<SetTriangle>    links_;
        Vertices *                  all_vertices_;
        VectPVertex                 vertices_;
        bool                        outermost_;
    };

    typedef std::vector<Mesh>        Meshes;
}

#endif  //  ! OPENMEEG_MESH_H
