// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include "vector.h"
#include "sparse_matrix.h"
#include "symmatrix.h"
#include "matrix.h"
#include "commandline.h"

#include <cmath>

using namespace OpenMEEG;

template <typename MATRIX>
void print_infos(const std::string& filename) {
    MATRIX M;
    M.load(filename);
    M.info();
}

int
main(int argc,char* argv[]) {

    print_version(argv[0]);

    //TODO doesn't say txt, if you don't specify it

    const CommandLine cmd(argc,argv,"Provides informations on a Matrix generated with OpenMEEG");

    const std::string& filename = cmd.option("-i",     std::string(),"Matrix file");
    const bool         sym      = cmd.option("-sym",   false,        "Data are symmetric matrices");
    const bool         sparse   = cmd.option("-sparse",false,        "Data are sparse matrices");
    
    if (cmd.help_mode())
        return 0;

    if (filename=="") {
        std::cerr << "Please set Matrix File" << std::endl;
        exit(1);
    }

    std::cout << "Loading : " << filename << std::endl;

    if (sym) {
        print_infos<SymMatrix>(filename);
    } else if (sparse) {
        print_infos<SparseMatrix>(filename);
    } else {
        print_infos<Matrix>(filename);
    }

    return 0;
}
