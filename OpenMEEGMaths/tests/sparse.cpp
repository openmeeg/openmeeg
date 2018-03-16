/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

The OpenMEEG software is a C++ package for solving the forward/inverse
problems of electroencephalography and magnetoencephalography.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's authors,  the holders of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

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
