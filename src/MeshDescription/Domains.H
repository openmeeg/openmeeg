/* FILE: $Id$ */

/*
Project Name : OpenMEEG

version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

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

#ifndef IMPLICITMESH_DOMAINS_H
#define IMPLICITMESH_DOMAINS_H

namespace MeshDescription {

    //  A domain is the association of a name and a vector of pairs.
    //  The first element of each pair corresponds to an interface.
    //  The second element of each pair states whether the domain is contains in the
    //  inside or the ouside of the volume bounded by the interface.

    typedef enum { Inside, Outside } InOut;
    typedef unsigned InterfaceId;

    //  A simple domain (HalfSpace) is given by an id (of type InterfaceId) identifying a closed surface and a side (of type InOut) information.
    //  The closed surface split the space into two components. The side depicts which of these two components is the simple domain.

    struct HalfSpace: protected std::pair<InterfaceId,InOut> {
        typedef std::pair<InterfaceId,InOut> base;
        HalfSpace(const int num): base(abs(num)-1,((num<0) ? Inside : Outside)) { }
        InterfaceId interface() const { return base::first;  }
        InOut       inout()     const { return base::second; }
    };

    //  A Domain is the intersection of simple domains (of type HalfSpace).
    //  In addition the domain is named.

    struct Domain: public std::vector<HalfSpace> {

        Domain() { }

        //  The name of the domain.

              std::string& name()       { return name_; }
        const std::string& name() const { return name_; }

    private:

        std::string name_;       // Name of the domain.
    };

    //  Domains is just a collection of Domain (here a simple vector).

    struct Domains: public std::vector<Domain> {
        unsigned index(const Domain& dom) const { return &dom-&*begin(); }
    };
}

#endif  //  ! IMPLICITMESH_DOMAINS_H
