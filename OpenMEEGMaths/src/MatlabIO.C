// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <MatlabIO.H>

namespace OpenMEEG {
    namespace maths {
        namespace details {
            const char helper<Vector>::message[] = "double Vector";
            const char helper<Matrix>::message[] = "2D full double Matrix";
            matvar_t* helper<SymMatrix>::saved;
            const char helper<SymMatrix>::message[] = "symmetric double Matrix";
            const char helper<SparseMatrix>::message[] = "2D full double Matrix";
        }
        const MatlabIO           MatlabIO::prototype;
        const std::string        MatlabIO::MagicTag("MATLAB");
        const MatlabIO::Suffixes MatlabIO::suffs = MatlabIO::init();
        const std::string        MatlabIO::Identity("matlab");
    }
}
