// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <cstring>

#include <matrix.h>
#include <symmatrix.h>
#include <vector.h>

#include <commandline.h>
#include <om_utils.h>

using namespace OpenMEEG;

void
help(const char* cmd_name) {
    std::cout << cmd_name <<" [-option] [filepaths...]" << std::endl << std::endl
              << "   Inverse HeadMatrix " << std::endl
              << "   Filepaths are in order :" << std::endl
              << "       HeadMat (bin), HeadMatInv (bin)" << std::endl << std::endl;

    exit(0);
}

int
main(int argc,char* argv[]) try {

    print_version(argv[0]);
    const CommandLine cmd(argc,argv);

    if (cmd.help_mode()) {
        help(argv[0]);
        return 0;
    }

    if (argc<3) {
        std::cerr << "Not enough arguments." << std::endl;
        help(argv[0]);
        return 1;
    }

    cmd.print();

    // Start Chrono

    const auto start_time = std::chrono::system_clock::now();

    SymMatrix HeadMat;

    HeadMat.load(argv[1]);
    HeadMat.invert(); // invert inplace
    HeadMat.save(argv[2]);

    // Stop Chrono

    const auto end_time = std::chrono::system_clock::now();
    dispEllapsed(end_time-start_time);

    return 0;
} catch (const OpenMEEG::maths::Exception& e) {
    std::cerr << e.what() << std::endl;
    return e.code();
} catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
    return 1;
}
