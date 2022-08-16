// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include "Properties.H"
#include "DataTag.H"

namespace OpenMEEG {

    template <typename REP=double>
    class OPENMEEG_EXPORT Conductivity {
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
}

namespace Types {
    template<>
    struct DataTrait<Utils::Properties::Named<std::string,OpenMEEG::Conductivity<double> > >{
        static const char TAG[];
    };
    const char DataTrait<Utils::Properties::Named<std::string,OpenMEEG::Conductivity<double> > >::TAG[]= "Conductivities";
}
