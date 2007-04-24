#include <iostream>

#include "mesh3.h"

int main (int argc, char** argv)
{
    if(argc==1) {
        std::cerr << "Wrong nb of parameters" << std::endl;
        exit(1);
    }
    mesh geo;
    geo.load(argv[1]);
    
    std::cerr << geo.nbr_pts() << ' ' << geo.nbr_trg() << std::endl;

    return 0;
}
