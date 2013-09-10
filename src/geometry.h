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
#include <vertex.h>
#include <triangle.h>
#include <mesh.h>
#include <interface.h>
#include <domain.h>

#include <vector>
#include <string>

namespace OpenMEEG {

    /** \brief Geometry contains the electrophysiological model
        Here are stored the vertices, meshes and domains
     */

    class OPENMEEG_EXPORT Geometry 
    {
        typedef enum { IDENTITY, INVERSE, INDICATOR} Function;

        /// friend class for reading geom/cond files.
        friend class GeometryReader;

    public:

        /// Default iterator of a Geometry is an Iterator on the meshes
        typedef Meshes::iterator          iterator;
        typedef Meshes::const_iterator    const_iterator;

        /// Iterators
        iterator                   begin()                 { return meshes_.begin();   }
        const_iterator             begin()           const { return meshes_.begin();   }
        iterator                   end()                   { return meshes_.end();     }
        const_iterator             end()             const { return meshes_.end();     }
        Vertices::iterator         vertex_begin()          { return vertices_.begin(); }
        Vertices::const_iterator   vertex_begin()    const { return vertices_.begin(); }
        Vertices::iterator         vertex_end()            { return vertices_.end();   }
        Vertices::const_iterator   vertex_end()      const { return vertices_.end();   }
        Domains::iterator          domain_begin()          { return domains_.begin();  }
        Domains::const_iterator    domain_begin()    const { return domains_.begin();  }
        Domains::iterator          domain_end()            { return domains_.end();    }
        Domains::const_iterator    domain_end()      const { return domains_.end();    }

        /// Constructors
        Geometry(): has_cond_(false), is_nested_(false), size_(0)  { }

              void       info()                           const; ///< \brief Print information on the geometry
        const bool&      has_cond()                       const { return has_cond_; }
        const bool&      is_nested()                      const { return is_nested_; }
        const bool       selfCheck()                      const; ///< \brief the geometry meshes intersect each other
        const bool       check(const Mesh& m)             const; ///< \brief check if m intersect geometry meshes
        const Vertices&  vertices()                       const { return vertices_; } ///< \brief returns the geometry vertices
        const Meshes&    meshes()                         const { return meshes_; } ///< \brief returns the geometry meshes
        const Domains&   domains()                        const { return domains_; } ///< \brief returns the geometry domains
        const unsigned   size()                           const { return size_; } ///< \brief the total number of vertices + triangles
        const unsigned   nb_vertices()                    const { return vertices_.size(); }
        const unsigned   nb_domains()                     const { return domains_.size(); }
        const unsigned   nb_meshes()                      const { return meshes_.size(); }
        const Interface& outermost_interface()            const; ///< \brief returns the outermost Interface of the geometry
        const Interface& interface(const std::string& id) const; ///< \brief returns the Interface called id \param id Interface name
        const Domain&    domain(const std::string&)       const; ///< \brief returns the Domain called id \param id Domain name
        const Domain&    domain(const Vect3& p)           const; ///< \brief returns the Domain containing the point p \param p a point

        void import_meshes(const Meshes& m); ///< \brief imports meshes from a list of meshes

        const double  sigma     (const Domain& d)                const { return (d.sigma()); }
        const double  sigma     (const Mesh& m1, const Mesh& m2) const { return funct_on_domains(m1, m2, IDENTITY); } // return the (sum) conductivity(ies) of the shared domain(s).
        const double  sigma_inv (const Mesh& m1, const Mesh& m2) const { return funct_on_domains(m1, m2, INVERSE); } // return the (sum) inverse of conductivity(ies) of the shared domain(s).
        const double  indicator (const Mesh& m1, const Mesh& m2) const { return funct_on_domains(m1, m2, INDICATOR); } // return the (sum) indicator function of the shared domain(s).
        const double  sigma_diff(const Mesh& m) const; // return the difference of conductivities of the 2 domains.
        const double  sigma     (const std::string&) const;
        const double  oriented(const Mesh&, const Mesh&) const;

        void read(const std::string& geomFileName, const std::string& condFileName = "", const bool OLD_ORDERING = false);
        void load_vtp(const std::string& filename);
        void write_vtp(const std::string& filename) const;

    private:

        Mesh& mesh(const std::string& id); ///< \brief returns the Mesh called id \param id Mesh name

        /// Members
        Vertices   vertices_;
        Meshes     meshes_;
        Domains    domains_;
        bool       has_cond_;
        bool       is_nested_;
        unsigned   size_;   // total number = nb of vertices + nb of triangles

        void          generate_indices(const bool);
        const Domains common_domains(const Mesh&, const Mesh&) const;
        const double  funct_on_domains(const Mesh&, const Mesh&, const Function& ) const;
    };
}

#endif  //! OPENMEEG_GEOMETRY_H
