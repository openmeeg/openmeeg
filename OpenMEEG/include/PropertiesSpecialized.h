// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include "Properties.H"
#include "DataTag.H"

#include <OpenMEEG_Export.h>

namespace OpenMEEG {

    class OPENMEEG_EXPORT ConductivityProperty {
    public:

        ConductivityProperty(): conductivity(1.0) { }

        double& value()       { return conductivity; }
        double  value() const { return conductivity; }

    private:

        double conductivity;    //  The conductivity of the layer (constant per layer).
    };

    inline std::istream& operator>>(std::istream& is,ConductivityProperty& m) { return is >> m.value(); }

    inline std::ostream& operator<<(std::ostream& os,const ConductivityProperty& m) { return os << m.value(); }
}

namespace Types {
    template<>
    struct DataTrait<Utils::Properties::Named<std::string,OpenMEEG::ConductivityProperty>>{
        static const char TAG[];
    };
    const char DataTrait<Utils::Properties::Named<std::string,OpenMEEG::ConductivityProperty>>::TAG[]= "Conductivities";
}
