// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <cmath>
#include <iostream>

#include <OpenMEEGMathsConfig.h>
#include <matrix.h>
#include <sparse_matrix.h>
#include <generic_test.hpp>

int main () {

    using namespace OpenMEEG;

    // section Matrix

    std::cout << std::endl << "========== matrices ==========" << std::endl;

    Matrix M(5,5);

    for (size_t i=0;i<M.nlin();++i)
        for (size_t j=0;j<M.ncol();++j)
            M(i,j) = pow(2.0,(double)i)+pow(2.0,(double)j);

    genericTest("full",M);

    Matrix Q = M.submat(3,1,2,3); // select submatrix
    std::cout << "Matrice Q : " << std::endl;
    Q.info();

    Matrix M1 = M;
    M1.insertmat(3,1,Q); // insert submatrix
    if (std::abs((M1-M).frobenius_norm()) > 1e-10) {
        std::cerr << "Error: insert matrix is WRONG" << std::endl;
        exit(1);
    }

    Matrix P(3,3);
    P(0,0) = 25 ; P(0,1) = 3 ; P(0,2) = 6;
    P(1,0) = 12 ; P(1,1) = 5 ; P(1,2) = 32;
    P(2,0) = 4  ; P(2,1) = 10; P(2,2) = 4;
    std::cout << "Matrice P : " << std::endl;
    P.info();

    Matrix Pinv = P.inverse();
    std::cout << "P Inverse Matrix : " << std::endl;
    Pinv.info();

    Matrix unit = P*Pinv;
    for (unsigned i=0;i<unit.nlin();++i) {
        for (unsigned j=0;j<unit.ncol();++j) {
            if (i==j) {
                if (std::abs(unit(i,j)-1)>eps) {
                    std::cerr << "Error: inverse is WRONG-1" << std::endl;
                    exit(1);
                }
            } else {
                if (std::abs(unit(i,j))>eps){
                    std::cerr << "Error: inverse is WRONG-2 " << "unit(" << i << "," << j << ") = " << unit(i,j) << std::endl;
                    exit(1);
                }
            }
        }
    }
    std::cout << std::endl << "BRAINVISA :" << std::endl;
    M.save("tmp.tex");
    M.load("tmp.tex");
    M.info();

    // SVD (wikipedia example)
    M1 = Matrix(4,5); M1.set(0.0);
    M1(0, 0) = 1; M1(0, 4) = 2;
    M1(1, 2) = 3; M1(3, 1) = 4;

    Matrix U, W;
    SparseMatrix S;
    M1.svd(U, S, W);
    std::cout << "SVD: M1 = U * S * W' " << std::endl;
    Matrix zero = M1 - U*S*W;
    if (zero.frobenius_norm()>eps) {
        std::cout << "M1 :" << std::endl;
        M1.info();
        std::cout << "U :" << std::endl;
        U.info();
        std::cout << "S :" << std::endl;
        S.info();
        std::cout << "W :" << std::endl;
        W.info();
        std::cerr << "Error: SVD is WRONG-1" << std::endl;
        exit(1);
    }

    // PseudoInverse
    M1.set(0.0);
    M1(0, 0) = 1; M1(0, 4) = 2;
    M1(1, 2) = 3; M1(3, 1) = 4;
    Matrix M1pinv = M1.pinverse();
    zero.set(0.);
    zero = M1*M1pinv*M1-M1;
    if (zero.frobenius_norm()>eps) {
        std::cout << "pinverse = M^{+} = " << std::endl;
        M1pinv.info();
        std::cout << "M*M^{+}*M-M = zero matrix = " << std::endl;
        zero.info();
        std::cerr << "Error: PseudoInverse is WRONG-2" << std::endl;
        exit(1);
    }
    return 0;
}
