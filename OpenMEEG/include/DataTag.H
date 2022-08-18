// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <string>
#include <iostream>
#include <IOUtils.H>
#include <OMExceptions.H>

namespace Types {

    //  This is just to be specialized. NEVER, EVER ATTEMPT TO INSTANCIATE IT.
    //  The actual class must containt a static string named TAG
    //  containing a tag identifying the type of data.

    template <typename T>
    struct DataTrait { };

    template <typename T>
    struct DataTag {
        static const char* tag() { return DataTrait<T>::TAG; }
    };

    //  The standard method just verifies that the tag is indeed correct.
    //  If not it raises an exception.
    //  This methods can be overloaded for more specific behaviours.

    template <typename T>
    std::istream& operator>>(std::istream& is,Types::DataTag<T>& tag) {
        bool DataP;
        is >> io_utils::match_optional(tag.tag(),DataP);

        if (!DataP)
            throw OpenMEEG::BadContent("data", tag.tag());

        return is;
    }
}
