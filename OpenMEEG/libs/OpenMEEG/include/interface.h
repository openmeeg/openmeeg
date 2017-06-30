/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

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

#pragma once

#include <string>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <om_common.h>
#include <mesh.h>

namespace OpenMEEG {

    /// An Oriented Mesh is a mesh associated with a boolean stating if it is well oriented. 
    class OrientedMesh: public std::pair<Mesh *, bool> 
    {
        typedef std::pair<Mesh *, bool> base;

    public:

        ///default constructors
        OrientedMesh() {}

        OrientedMesh(Mesh& _mesh, bool _inside): base(&_mesh, _inside) {}

              Mesh&  mesh()              { return *first;  } ///< \brief access mesh
        const Mesh&  mesh()        const { return *first;  } ///< \brief access mesh
              int    orientation() const { return ( second )?1:-1; } ///< \brief orientation is +1 or -1 ?
    };

    /** Interface class
        An interface is a closed-shape composed of oriented meshes (here pointer to meshes)
    **/
    class OPENMEEG_EXPORT Interface: public std::vector<OrientedMesh> {

    public:

        typedef Mesh::VectPTriangle VectPTriangle;

        /// Default Constructor
        Interface(): name_(""), outermost_(false) { }
        
        /// Constructor from a name
        Interface(const std::string _name): name_(_name), outermost_(false) { }

        const std::string   name()                       const      { return name_; } ///< \return Interface name
        const bool &        outermost()                  const      { return outermost_; } ///< \return true if it is the outermost interface.
              void          set_to_outermost(); ///< set all interface meshes to outermost state.
              bool          contains_point(const Vect3& p) const; ///< \param p a point \return true if point is inside interface
              bool          check(bool checked = false); ///< Check the global orientation

        /// \return the total number of the interface vertices
        size_t nb_vertices() const {
            size_t nb = 0;
            for ( const_iterator omit = begin(); omit != end(); ++omit) {
                nb += omit->mesh().nb_vertices();
            }
            return nb;
        }
        
        /// \return the total number of the interface triangles
        size_t nb_triangles() const {
            size_t nb = 0;
            for ( const_iterator omit = begin(); omit != end(); ++omit) {
                nb += omit->mesh().nb_triangles();
            }
            return nb;
        }

        /// \return the adjacent triangles
        VectPTriangle adjacent_triangles(const Triangle& t) const {
            VectPTriangle triangles;
            for ( const_iterator omit = begin(); omit != end(); ++omit) {
                VectPTriangle tri = omit->mesh().adjacent_triangles(t);
                triangles.insert(triangles.end(), tri.begin(), tri.end());
            }
            return triangles;
        }

    private:

        double compute_solid_angle(const Vect3& p) const; ///< Given a point p, it computes the solid angle \return should return +/- 4 PI or 0.

        std::string name_;      ///< is "" by default
        bool        outermost_; ///< tell weather or not the interface touches the Air (Outermost) Domain.
    };

    /// A vector of Interface is called Interfaces
    typedef std::vector<Interface> Interfaces;
}
