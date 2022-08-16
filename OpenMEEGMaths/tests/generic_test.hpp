// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <iostream>

#include <vector.h>

double eps = 1.e-12;

template <typename T>
void genericTest(const std::string& basename, T &M) {

    using namespace OpenMEEG;

    std::cout << " Generic Test " << std::endl;
    std::cout << "   nlin  = " << static_cast<int>(M.nlin()) << std::endl;
    std::cout << "   ncol  = " << static_cast<int>(M.ncol()) << std::endl;
    Vector v(M.ncol());
    v.set(1);
    v = M*v;

    std::cout << std::endl << "BASE :" << std::endl;
    M.info();

    // Test IO
    std::cout << std::endl << "BIN :" << std::endl;
    const std::string binname = basename+".bin";
    M.save(binname);
    M.load(binname);
    M.info();

    std::cout << std::endl << "TXT :" << std::endl;
    const std::string txtname = basename+".txt";
    M.save(txtname);
    M.load(txtname);
    M.info();

    // TODO Here for sparse matrix 

    const std::string matname = basename+".mat";
    std::cout << "MAT :" << std::endl;
    M.save(matname);
    M.load(matname);
    M.info();

    std::cout << "   operator * OK" << std::endl;
    std::cout.flush();
}
