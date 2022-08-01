/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

The OpenMEEG software is a C++ package for solving the forward/inverse
problems of electroencephalography and magnetoencephalography.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's authors,  the holders of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

#include <om_utils.h>
#include <commandline.h>
#include <forward.h>

using namespace OpenMEEG;

void
getHelp(const char* command) {
    std::cout << command << " [-h | --help] filepaths" << std::endl << std::endl
              << "   Compute the forward problem " << std::endl
              << "   Filepaths are in order :" << std::endl
              << "   GainMatrix (bin), RealSourcesData (txt), SimulatedData (txt), NoiseLevel (float)" << std::endl
              << std::endl;
}

void error(const char* command,const bool unknown_option=false) {
    std::cerr << "Error: " << ((unknown_option) ? "Unknown option." : "Not enough arguments.") << std::endl;
    getHelp(command);
    exit(1);
}

int
main(int argc,char **argv) {

    print_version(argv[0]);

    if (argc==2 && (!strcmp(argv[1],"-h") || !strcmp(argv[1],"--help"))) {
        getHelp(argv[0]);
        return 0;
    }

    if (argc<5)
        error(argv[0]);

    // Start Chrono

    auto start_time = std::chrono::system_clock::now();

    print_commandline(argc,argv);

    // declaration of argument variables======================================================================

    Matrix GainMatrix(argv[1]);
    Matrix RealSourcesData(argv[2]);

    const double NoiseLevel = atof(argv[4]);

    Forward SimulatedData(GainMatrix,RealSourcesData,NoiseLevel);

    // write output variables ===================================================================================

    SimulatedData.save(argv[3]);

    // Stop Chrono

    auto end_time = std::chrono::system_clock::now();
    dispEllapsed(end_time-start_time);

    return 0;
}
