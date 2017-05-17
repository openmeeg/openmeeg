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
