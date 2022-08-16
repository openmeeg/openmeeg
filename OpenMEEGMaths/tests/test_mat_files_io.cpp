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

int main (int argc, char** argv) {

    using namespace OpenMEEG;

    if ( argc != 3 ) 
    {
        std::cerr << "Wrong nb of parameters. Should be 2." << std::endl;
        exit(1);
    }

    std::cout << "Mesh : " << argv[1] << std::endl;

    try {
        Matrix matrix;
        matrix.load(argv[1]);
        matrix.info();

        matrix.load(argv[2]);
        matrix.info();
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        exit(1);
    }

}
