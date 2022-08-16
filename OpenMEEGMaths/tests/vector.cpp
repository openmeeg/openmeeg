// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <iostream>

#include <OpenMEEGMathsConfig.h>
#include <vector.h>
#include <matrix.h>
#include <generic_test.hpp>

int main () {

    using namespace OpenMEEG;

    // section Vector

    std::cout << std::endl << "========== vectors ==========" << std::endl;

    Vector v(8);
    v.set(0);
    v.save("tmp.bin");

    for (int i=0;i<8;++i)
        v(i) = i;

    v.save("tmp.txt");
    v.load("tmp.bin");
    std::cout << "v = " << std::endl << v << std::endl;
    v.load("tmp.txt");
    std::cout << "v = " << std::endl << v << std::endl;

    std::cout << "MAT :" << std::endl;
    v.save("tmp_matrix.mat");
    v.load("tmp_matrix.mat");
    v.info();

    v(0) = 115;
    v(7) = 0.16;
    v(3) = 0.22;
    v(2) = 2.;
    Matrix m(v,v.size(),1);
    m = m* m.transpose();
    m -= v.outer_product(v);
    if ( m.frobenius_norm() > eps) {
        m.info();
        std::cerr << "Error: Vector outerproduct is WRONG-1" << std::endl;
        exit(1);
    }

    return 0;
}
