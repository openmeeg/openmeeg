/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
*/

#ifndef UTILS_PROPERTIES_SPECIALIZED_H
#define UTILS_PROPERTIES_SPECIALIZED_H

#include "Properties.H"
#include "DataTag.H"

template <typename REP=double>
class Conductivity {
public:

    Conductivity(): conductivity(1.0) { }

    REP& sigma()       { return conductivity; }
    REP  sigma() const { return conductivity; }
private:
    REP conductivity;    //  The conductivity of the layer (constant per layer).
};

template <typename REP>
inline std::istream& operator>>(std::istream& is,Conductivity<REP>& m) { return is >> m.sigma(); }

template <typename REP>
inline std::ostream& operator<<(std::ostream& os,const Conductivity<REP>& m) { return os << m.sigma(); }

namespace Types {
    template<>
    struct DataTrait<Utils::Properties::Named<std::string, Conductivity<double> > >{
        static const char TAG[];
    };
    const char DataTrait<Utils::Properties::Named<std::string, Conductivity<double> > >::TAG[]= "Conductivities";
};


#endif  //  ! UTILS_PROPERTIES_SPECIALIZED_H
