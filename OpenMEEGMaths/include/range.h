// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <iostream>
#include <vector>

namespace OpenMEEG::maths {

    class Range {
    public:

        Range(): start_index(0),end_index(0) { }
        Range(const size_t s,const size_t e): start_index(s),end_index(e) { }

        size_t start() const { return start_index; }
        size_t end()   const { return end_index;   }

        size_t length() const { return end()-start()+1; }

        bool contains(const size_t ind) const { return ind>=start() && ind<=end();               }
        bool intersect(const Range& r)  const { return contains(r.start()) || contains(r.end()); }

        bool operator==(const Range& r) const { return start()==r.start() && end()==r.end();     }
        bool operator!=(const Range& r) const { return start()!=r.start() || end()!=r.end();     }

    private:

        size_t start_index;
        size_t end_index;
    };

    inline std::ostream& operator<<(std::ostream& os,const Range& r) {
        return os << '(' << r.start() << ',' << r.end() << ')';
    }
}
