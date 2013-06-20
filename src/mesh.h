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

#include <vector>
#include <set>
#include <string>
#include "triangle.h"

namespace OpenMEEG {

    /** \brief  Mesh

        Mesh class

    **/

    class OPENMEEG_EXPORT Mesh: public Triangles {

    public:
        
        typedef std::set<Vertex *, std::less<Vertex *> >    SetPVertex;
        typedef std::vector<Vertex *>                       VectPVertex;
        typedef VectPVertex::iterator                       vertex_iterator;
        typedef VectPVertex::const_iterator                 const_vertex_iterator;

        Mesh(): _name(""), _vertices(NULL), _outermost(false) {};

        Mesh(Vertices *vertices, const std::string name = ""): _vertices(vertices), _name(name), _outermost(false) { }

        Mesh(std::string name): _name(name), _vertices(NULL), _outermost(false) { }

        // Iterators on vertices
              vertex_iterator      vertex_begin()                { return mesh_vertices().begin(); }
        const_vertex_iterator      vertex_begin()  const         { return mesh_vertices().begin(); }
              vertex_iterator      vertex_end()                  { return mesh_vertices().end(); }
        const_vertex_iterator      vertex_end()    const         { return mesh_vertices().end(); }

        std::vector<SetTriangle>   links()         const         { return _links; }
        std::vector<SetTriangle> & links()                       { return _links; }

        std::string                name()          const         { return _name; }
        std::string &              name()                        { return _name; }

        VectPVertex                mesh_vertices() const         { return _mesh_vertices; }
        VectPVertex &              mesh_vertices()               { return _mesh_vertices; }

        const size_t               nb_vertices()   const         { return _mesh_vertices.size(); }
        const size_t               nb_triangles()  const         { return this->size(); }

        const Vertices *           vertices()      const         { return _vertices; }
              Vertices *           vertices()                    { return _vertices; }

        // mesh state
        void mesh_info() const ;
        bool triangle_intersection( const Mesh, const Triangle, const Mesh, const Triangle) const;
        bool has_self_intersection() const;
        bool intersection(const Mesh) const;
        bool has_correct_orientation() const;
        void update();

        const SetTriangle& get_triangles_for_point(const Vertex& V) const {
            size_t i = 0;
            for (const_vertex_iterator vit = mesh_vertices().begin(); vit != mesh_vertices().end(); vit++, i++) {
                if (*vit == &V) {
                    return links()[i];
                }
            }
        }

        //  Returns True if it is an outermost mesh.
              bool&        outermost()       { return _outermost; }
        const bool&        outermost() const { return _outermost; }

    private:
        
        std::string                 _name;
        std::vector<SetTriangle>    _links;
        Vertices *                  _vertices;
        VectPVertex                 _mesh_vertices;
        bool                        _outermost;

    };

    typedef std::vector<Mesh>        Meshes;
}

#endif  //  ! OPENMEEG_MESH_H
