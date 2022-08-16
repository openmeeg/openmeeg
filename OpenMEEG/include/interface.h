// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <om_common.h>
#include <mesh.h>

namespace OpenMEEG {

    class Mesh;

    /// An Oriented Mesh is a mesh associated with a boolean stating if it is well oriented.

    class OrientedMesh {
    public:

        typedef enum { Normal=1, Opposite=-1 } Orientation ;

        OrientedMesh() {}

        OrientedMesh(Mesh& m,const Orientation o): meshptr(&m),orient(o) { }

              Mesh&  mesh()       { return *meshptr; } ///< \brief access mesh
        const Mesh&  mesh() const { return *meshptr; } ///< \brief access mesh

        int orientation() const { return orient; } ///< \brief orientation is +1 or -1 ?

        void change_orientation() { orient = -orient; }

    private:

        Mesh* meshptr;
        int   orient;
    };

    /// Interface class
    /// An interface is a closed-shape composed of oriented meshes (vector of oriented meshes).

    class OPENMEEG_EXPORT Interface {
    public:

        //        using OrientedMesh::Orientation;

        typedef std::vector<OrientedMesh> OrientedMeshes;

        /// Default Constructor

        Interface() { }

        /// Constructor from a name

        Interface(const std::string& interfname): interface_name(interfname) { }

        const std::string& name() const { return interface_name; } ///< \return Interface name

              OrientedMeshes& oriented_meshes()       { return orientedmeshes; }
        const OrientedMeshes& oriented_meshes() const { return orientedmeshes; }

        bool outermost() const { return outermost_interface; } ///< \return true if it is the outermost interface.
        void set_to_outermost(); ///< set all interface meshes to outermost state.

        bool contains(const Vect3& p) const; ///< \param p a point \return true if point is inside interface

        bool is_mesh_orientations_coherent(const bool doublechecked=false); ///< Check the global orientation

        /// \return the total number of the interface vertices

        size_t nb_vertices() const {
            size_t nb = 0;
            for (const auto& omesh : oriented_meshes())
                nb += omesh.mesh().vertices().size();
            return nb;
        }

        /// \return the total number of the interface triangles

        size_t nb_triangles() const {
            size_t nb = 0;
            for (const auto& omesh : oriented_meshes())
                nb += omesh.mesh().triangles().size();
            return nb;
        }

        /// \return the adjacent triangles

        TrianglesRefs adjacent_triangles(const Triangle& t) const {
            TrianglesRefs triangles;
            for (const auto& omesh : oriented_meshes()) {
                const TrianglesRefs& tri = omesh.mesh().adjacent_triangles(t);
                triangles.insert(triangles.end(),tri.begin(),tri.end());
            }
            return triangles;
        }

    private:

        double solid_angle(const Vect3& p) const; ///< Given a point p, it computes the solid angle \return should return +/- 4 PI or 0.

        std::string    interface_name      = "";    ///< interface name is "" by default
        bool           outermost_interface = false; ///< whether or not the interface touches the Air (outermost) domain.
        OrientedMeshes orientedmeshes;
    };

    /// A vector of Interface is called Interfaces.

    typedef std::vector<Interface> Interfaces;
}
