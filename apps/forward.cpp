// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <om_utils.h>
#include <commandline.h>
#include <forward.h>

using namespace OpenMEEG;

void
help(const char* command) {
    std::cout << command << " [-h | --help] filepaths" << std::endl << std::endl
              << "   Compute the forward problem " << std::endl
              << "   Filepaths are in order :" << std::endl
              << "   GainMatrix (bin), RealSourcesData (txt), SimulatedData (txt), NoiseLevel (float)" << std::endl
              << std::endl;
}

int
main(int argc,char **argv) {

    print_version(argv[0]);
    const CommandLine cmd(argc,argv);

    if (cmd.help_mode()) {
        help(argv[0]);
        return 0;
    }

    if (argc<5) {
        std::cerr << "Not enough arguments." << std::endl;
        help(argv[0]);
        return 1;
    }

    print_version(argv[0]);
    cmd.print();

    // Start Chrono

    const auto start_time = std::chrono::system_clock::now();

    // Loading input matrices.

    const Matrix GainMatrix(argv[1]);
    const Matrix RealSourcesData(argv[2]);

    const double NoiseLevel = atof(argv[4]);

    const Forward SimulatedData(GainMatrix,RealSourcesData,NoiseLevel);

    // Write output variables

    SimulatedData.save(argv[3]);

    // Stop Chrono

    const auto end_time = std::chrono::system_clock::now();
    dispEllapsed(end_time-start_time);

    return 0;
}
