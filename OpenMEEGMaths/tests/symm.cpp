// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <cmath>
#include <iostream>

#include <OpenMEEGMathsConfig.h>
#include <symmatrix.h>
#include <matrix.h>
#include <generic_test.hpp>

int main() {

    using namespace OpenMEEG;

    // section SymMatrix

    std::cout << std::endl << "========== symmetric matrices ==========" << std::endl;

    SymMatrix S(4);
    for (unsigned i=0;i<4;++i)
        for (unsigned j=i; j<4;++j)
            S(i,j) = pow(2.0,(double)i)+pow(3.0,(double)j);

    genericTest("symm",S);

    const Matrix R = S(1,2,0,2); // extract submatrix
    std::cout << "Matrice R : " << std::endl;
    R.info();

    return 0;
}
