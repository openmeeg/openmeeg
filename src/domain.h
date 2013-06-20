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

#ifndef OPENMEEG_DOMAIN_H
#define OPENMEEG_DOMAIN_H

#include <string>
#include "interface.h"

namespace OpenMEEG {

    //  A domain is the association of a name and a vector of pairs.
    //  The first element of each pair corresponds to an interface.
    //  The second element of each pair states whether the domain is contains in the
    //  inside or the ouside of the volume bounded by the interface.

    typedef enum { Inside, Outside } InOut;

    //  A simple domain (HalfSpace) is given by an interface (of type Interface) identifying a closed surface and a side (of type InOut) information.
    //  The closed surface split the space into two components. The side depicts which of these two components is the simple domain.

    class HalfSpace: private std::pair<const Interface *, bool> {

        typedef std::pair<const Interface *, bool> base;

    public:

        HalfSpace(const Interface * interface, bool inside) { base(std::make_pair<const Interface *, bool>(interface, inside)); }

        Interface interface() const { return (*base::first);  }
        bool      inside()    const { return (base::second); }
    };

    //  A Domain is the intersection of simple domains (of type HalfSpace).
    //  In addition the domain is named, has conductivity

    class Domain: public std::vector<HalfSpace> {

        typedef std::vector<HalfSpace> base;

    public:

        Domain(): _name(""), _conductivity(0.), _innermost(false) { }

        //  The interfaces of the domain.
        // const std::vector<Interface *> interfaces() const  {
            // std::vector<Interface *> _interfaces(this->size());
            // size_t i = 0;
            // for (base::const_iterator hit = this->begin(); hit != this->end(); hit++, i++) {
                // _interfaces[i] = this->first;
            // }
            // return _interfaces;
        // }
        
        //  The name of the domain.
              std::string& name()            { return _name; }
        const std::string& name()      const { return _name; }
        
        //  The conductivity of the domain.
              double&      sigma()           { return _conductivity; }
        const double&      sigma()     const { return _conductivity; }

        //  Returns the innermost state of the domain.
              bool&        innermost()       { return _innermost; }
        const bool&        innermost() const { return _innermost; }

        bool contains_point(const Vect3&) const;

        int meshOrient(const Mesh& m) { 
            // return 1 if the mesh is oriented toward the domain
                  // -1 if not
                  //  0 else (the mesh is not part of the domain boundary)
            // for (Interface::const_iterator mit = interface().begin(); mit != interface().end(); mit++) {
                // if (mit == m ) {
                    // return (mit.second)?1:-1;
                // }
            // }
            return 0;
        }
    private:

        std::string _name;         // Name of the domain.
        double      _conductivity; // Conductivity of the domain.
        bool        _innermost;    // Innermost domain ?
    };

    bool operator==(const Domain& d1, const Domain& d2) {
        for (Domain::const_iterator dit1 = d1.begin(), dit2 = d2.begin(); dit1 != d1.end(); ++dit1, ++dit2) {
            if (&*dit1 != &*dit2) {
                return false;
            }
        }
        return true;
    }

    typedef std::vector<Domain >     Domains;

    /*
        class OPENMEEG_EXPORT Domain {
            // A domain contains information about its surrounding meshes (which must defines a closed-shape) and its conductivity sigma
            public:
                Domain(char * _name, Map_MeshOriented _map_mesh, double _conductivity = 0.)    {
                    this->nam()          = _name;
                    this->conductivity() = _conductivity;
                    this->map_mesh       = _map_mesh;
                }
                const   double  sigma()         {return conductivity;}
                        char *  name()          {return nam;}
                        int     meshOrient(Mesh * m);

                // bool operator< (Domain d) {return (strcmp(this-> name(), d.name()) < 0); }

            private:
                double           conductivity;
                char *           nam;
                Map_MeshOriented map_mesh;

                int meshOrient(Mesh * m) { // return 1 if the mesh is oriented toward the domain
                                                 // -1 if not
                                                 //  0 else (the mesh is not part of the domain boundary)
                    Map_MeshOriented::const_iterator mit = map_mesh.find(m);
                    if (mit != map_mesh.end() ) {
                        return (mit.second)?1:-1;
                    }
                    return 0;
                }
        };
        */
}

#endif  //  ! OPENMEEG_DOMAIN_H
