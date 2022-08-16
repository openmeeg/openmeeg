// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <fstream>
#include <cstring>

#include <sensors.h>

using namespace OpenMEEG;

int main(const int,const char** argv) {
// usage : sensors sensors_file_description.txt

    /*** tests on sensors file ****/
    Sensors S(argv[1]);

    size_t n = S.getNumberOfSensors();
    std::cout << "Number of sensors of S : " << n << std::endl;

    if (S.isEmpty())
        std::cout << "WARNING : empty sensors !" << std::endl;
    else{
        S.info();

        /**** test on copy constructor ****/
        Sensors Scopy(S);
        if(Scopy.getNumberOfSensors()!=n)
            std::cout << "ERROR in copy from copy constructor : incorrect number of sensors" << std::endl;

        Scopy.info();
    }
}
