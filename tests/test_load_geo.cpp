#include <iostream>

#include "geometry.h"

using namespace OpenMEEG;

int
main(int argc,char** argv) {

	if (argc!=3) {
        std::cerr << "Wrong nb of parameters" << std::endl;
        return 1;
    }

    std::cerr << "Geom : " << argv[1] << std::endl
              << "Cond : " << argv[2] << std::endl;

    Geometry geo(argv[1],argv[2]);

    std::cerr << "Geometry degrees of freedom: " << geo.nb_parameters() << std::endl;

    return 0;
}
