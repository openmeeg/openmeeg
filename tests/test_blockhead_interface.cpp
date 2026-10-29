#include <iostream>

#include <PropertiesSpecialized.h>

#include <geometry.h>
#include <assemble.h>

using namespace OpenMEEG;

int
main(int argc,char** argv) {

    std::cerr << argc << std::endl;
	if (argc!=4) {
        std::cerr << "Wrong nb of parameters" << std::endl << std::endl
                  << "Usage: " << argv[0] << " geom_filename cond_filename output_filename" << std::endl;
        return 1;
    }

    std::cerr << "Geom : " << argv[1] << std::endl
              << "Conductivities : " << argv[2] << std::endl;

    Geometry geom(argv[1],argv[2]);
    geom.set_conductivity("Scalp",-1.0);
    geom.set_conductivity("Skull",-1.0);

    std::cerr << "Geometry degrees of freedom: " << geom.num_unknowns() << std::endl;

    const HeadMatrix& block_head_matrix = BlockHeadMat(geom);

    geom.set_conductivity("Scalp",1.0);
    geom.set_conductivity("Skull",0.0125);

    SymMatrix headmat = block_head_matrix.make_head_matrix();
    HeadMatrix::deflate(geom,headmat);
    headmat.save(argv[3]);

    return 0;
}
