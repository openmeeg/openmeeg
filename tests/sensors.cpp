#include <fstream>
#include <cstring>

#include <matrice.h>
#include <sensors.h>

int main(int argc, char** argv) {
// usage : sensors sensors_file_description.txt

    /*** tests on sensors file ****/
    Sensors S(argv[1]);

    size_t n = S.getNumberOfSensors();
    std::cout << "Number of sensors of S : " << n << std::endl;

    if (S.isEmpty())
        std::cout << "WARNING : empty sensors !" << std::endl;
    else{
        matrice& p = S.getPositions();
        matrice& o = S.getOrientations();
        std::cout << std::endl << "Positions of sensors : " << std::endl;
        std::cout << p << std::endl;
        std::cout << std::endl << "Orientations of sensors : " << std::endl;
        std::cout << o << std::endl;

        /**** test on copy constructor ****/
        Sensors Scopy(S);
        if(Scopy.getNumberOfSensors()!=n)
            std::cout << "ERROR in copy from copy constructor : incorrect number of sensors" << std::endl;

        std::cout << Scopy << std::endl;
    }
}
