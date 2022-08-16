// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <AsciiIO.H>

namespace OpenMEEG {
    namespace maths {
        const AsciiIO           AsciiIO::prototype;
        const AsciiIO::Suffixes AsciiIO::suffs = AsciiIO::init();
        const std::string       AsciiIO::Identity("ascii");
    }
}
