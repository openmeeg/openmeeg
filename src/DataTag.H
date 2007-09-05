#ifndef DATATAG_H
#define DATATAG_H

#include <string>
#include <iostream>
#include <IOUtils.H>
#include <FileExceptions.H>

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

#if WIN32
#pragma warning (disable: 4290)	// this particular throw is not the one Visual 2005 is waiting for
#endif
    template <typename T>
    std::istream& operator>>(std::istream& is,Types::DataTag<T>& tag) throw(std::io_except::read_error) {
#if WIN32
#pragma warning (default: 4290)
#endif
        bool DataP;
        is >> io_utils::match_optional(tag.tag(),DataP);

        if (!DataP)
            throw std::io_except::read_error(is);

        return is;
    }
}

#endif  // !DATATAG_H
