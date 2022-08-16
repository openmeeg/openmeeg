// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <iostream>

#include <OpenMEEGMathsConfig.h>
#include <sparse_matrix.h>
#include <fast_sparse_matrix.h>
#include <generic_test.hpp>

int main () {

    using namespace OpenMEEG;

    // section SparseMatrix

    std::cout << std::endl << "========== sparse matrices ==========" << std::endl;

    SparseMatrix spM(10,10);
    unsigned n = 0;
    for ( unsigned i=0;i<5;++i) {
        n = (n*1237+1493)%1723;
        const int p = (n*1237+1493)%1723;
        spM(n%10, p%10) = n;
    }
    genericTest("sparse",spM);

    Matrix U(10,10);
    U.set(1.0);
    U(2,3)=0.12;
    U(1,9)=12.01;
    U(4,8)=-2.1;
    Vector v(10);
    v.set(2.);
    v(3)=v(8)=0.11;
    // Mat & Sparse
    Matrix Mzero = spM*U - Matrix(spM)*U - Matrix(spM)*U + spM*U;
    // Sparse & Sparse
    SparseMatrix spM2(10,10);
    for ( unsigned i=0;i<5;++i) {
        n = (n*1007+1493)%2551;
        const int p = (n*1007+1493)%2551;
        spM2(n%10, p%10) = n;
    }
    Mzero += Matrix(spM*spM2) - Matrix(spM)*Matrix(spM2) - Matrix(spM2)*Matrix(spM) + Matrix(spM2*spM);
    // Vectt & Sparse
    Vector Vzero = (spM*v) - (Matrix(spM)*v);
    if ( Mzero.frobenius_norm() + Vzero.norm() > eps) {
        std::cerr << "Error: Operator* is WRONG-1" << std::endl;
        Mzero.info();
        Vzero.info();
        exit(1);
    }

    std::cout << std::endl << "========== fast sparse matrices ==========" << std::endl;
    std::cout << spM;
    FastSparseMatrix fspM(spM);
    std::cout << fspM;

    return 0;
}
