/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
lastrevision      : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
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

    //  A simple domain (SimpleDomain) is given by an id (of type InterfaceId) identifying a closed surface and a side (of type InOut) information.
    //  The closed surface split the space into two components. The side depicts which of these two components is the simple domain.

    struct SimpleDomain: protected std::pair<InterfaceId,InOut> {
        typedef std::pair<InterfaceId,InOut> base;
        SimpleDomain(const int num): base(abs(num)-1,((num<0) ? Inside : Outside)) { }
        InterfaceId interface() const { return base::first;  }
        InOut       inout()     const { return base::second; }
    };

    //  A Domain is the intersection of simple domains (of type SimpleDomain).
    //  In addition the domain is named.

    struct Domain: public std::vector<SimpleDomain> {

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
