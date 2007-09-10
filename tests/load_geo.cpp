#include <iostream>

#include "geometry.h"

int main (int argc, char** argv)
{
	if(argc!=3) {
        std::cerr << "Wrong nb of parameters" << std::endl;
        exit(1);
    }

    std::cerr << "Geom : " << argv[1] << std::endl;
    std::cerr << "Cond : " << argv[2] << std::endl;

    Geometry geo;
    int taille = geo.read(argv[1],argv[2]);

    std::cerr << "Geometry Size : " << taille << std::endl;

    return 0;
}
