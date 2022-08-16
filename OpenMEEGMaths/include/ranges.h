// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <iostream>
#include <vector>

#include <range.h>
#include <Exceptions.H>

#include "OpenMEEGMathsConfig.h"

namespace OpenMEEG::maths {

    class Ranges: public std::vector<Range> {

        typedef std::vector<Range> base;

    public:

        using base::base;

        unsigned add(const Range& range) {
            for (unsigned i=0;i<size();++i)
                if ((*this)[i].intersect(range)) {
                    if (range!=(*this)[i])
                        throw OverlappingRanges(range,(*this)[i]);
                    return i;
                }
            push_back(range);
            return size()-1;
        }

        unsigned find_index(const size_t ind) const {
            for (unsigned i=0;i<size();++i)
                if ((*this)[i].contains(ind))
                    return i;
            throw NonExistingBlock(ind);
        }

        unsigned find_index(const Range& range) const {
            for (unsigned i=0;i<size();++i)
                if ((*this)[i].intersect(range)) {
                    if (range!=(*this)[i])
                        throw OverlappingRanges(range,(*this)[i]);
                    return i;
                }
            throw NonExistingRange(range);
        }
    };
}
