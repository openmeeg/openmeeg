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
#include <gain.h>

using namespace OpenMEEG;

void getHelp(char** argv);

inline void
error(const char* command,const bool unknown_option=false) {
    std::cerr << "Error: " << ((unknown_option) ? "Unknown option." : "Not enough arguments.") << std::endl
              << "Please try \"" << command << " -h\" or \"" << command << " --help \" \n" << std::endl;
    exit(1);
}

int
main(int argc,char** argv) {

    print_version(argv[0]);

    if (argc<2)
        error(argv[0]);

    if ((!strcmp(argv[1],"-h")) || (!strcmp(argv[1],"--help")))
        getHelp(argv);

    // Start Chrono

    print_commandline(argc,argv);

    const std::string& option = argv[1];
    if (argc<5)
        error(argv[0]);

    const auto start_time = std::chrono::system_clock::now();

    if (!strcmp(argv[1],"-EEG")) {

        // EEG DATA

        if (argc<6)
            error(argv[0]);

        //  Split the 2 matrix multiplications in order to spare memory.
        //  This is why we do not use GainEEG...

        const SymMatrix HeadMatInv(argv[2]);
        const SparseMatrix Head2EEGMat(argv[4]);
        const Matrix& tmp = Head2EEGMat*HeadMatInv;
        const Matrix SourceMat(argv[3]);
        const Matrix& EEGGainMat = tmp*SourceMat;
        EEGGainMat.save(argv[5]);

    } else if (!strcmp(argv[1],"-EEGadjoint")) {

        // Compute gain matrix with the adjoint method for use with EEG DATA

        if (argc<8)
            error(argv[0]);

        Geometry geo(argv[2],argv[3]);
        const Matrix dipoles(argv[4]);
        const SymMatrix HeadMat(argv[5]);
        const SparseMatrix Head2EEGMat(argv[6]);

        const GainEEGadjoint EEGGainMat(geo, dipoles, HeadMat, Head2EEGMat);
        EEGGainMat.save(argv[7]);

    } else if (!strcmp(argv[1],"-MEG")) {

        // MEG DATA

        if (argc<7)
            error(argv[0]);

        //  We split the 3 matrix multiplications in order to spare memory.
        //  This is also why we do not use GainMEG...

        const SymMatrix HeadMatInv(argv[2]);
        const Matrix Head2MEGMat(argv[4]);
        const Matrix& tmp1 = Head2MEGMat*HeadMatInv;
        const Matrix SourceMat(argv[3]);
        const Matrix& tmp2 = tmp1*SourceMat;
        const Matrix Source2MEGMat(argv[5]);
        const Matrix MEGGainMat = Source2MEGMat+tmp2;
        MEGGainMat.save(argv[6]);

    } else if ( !strcmp(argv[1], "-MEGadjoint") ) {

        // Compute the gain matrix with the adjoint method for use with MEG DATA

        if (argc<9)
            error(argv[0]);

        Geometry geo(argv[2],argv[3]);
        const Matrix dipoles(argv[4]);
        const SymMatrix HeadMat(argv[5]);
        const Matrix Head2MEGMat(argv[6]);
        const Matrix Source2MEGMat(argv[7]);

        const GainMEGadjoint MEGGainMat(geo, dipoles, HeadMat, Head2MEGMat, Source2MEGMat);
        MEGGainMat.save(argv[8]);

    } else if (!strcmp(argv[1],"-EEGMEGadjoint")) {

        // Compute the gain matrices with the adjoint method for EEG and MEG DATA

        if (argc<11)
            error(argv[0]);

        Geometry geo(argv[2],argv[3]);
        const Matrix dipoles(argv[4]);
        const SymMatrix HeadMat(argv[5]);
        const SparseMatrix Head2EEGMat(argv[6]);
        const Matrix Head2MEGMat(argv[7]);
        const Matrix Source2MEGMat(argv[8]);

        const GainEEGMEGadjoint EEGMEGGainMat(geo,dipoles,HeadMat,Head2EEGMat,Head2MEGMat,Source2MEGMat);
        EEGMEGGainMat.saveEEG(argv[9]);
        EEGMEGGainMat.saveMEG(argv[10]);

    } else if (!strcmp(argv[1],"-InternalPotential") || !strcmp(argv[1],"-IP")) {

        if (argc<7)
            error(argv[0]);

        const SymMatrix HeadMatInv(argv[2]);
        const Matrix Head2IPMat(argv[4]);

        const Matrix& tmp1 = Head2IPMat*HeadMatInv;
        const Matrix SourceMat(argv[3]);
        const Matrix& tmp2 = tmp1*SourceMat;
        const Matrix Source2IPMat(argv[5]);

        const Matrix& InternalPotGainMat = Source2IPMat+tmp2;
        InternalPotGainMat.save(argv[6]);

    } else if (!strcmp(argv[1],"-EITInternalPotential") || !strcmp(argv[1], "-EITIP")) {

        if (argc<6)
            error(argv[0]);

        const SymMatrix HeadMatInv(argv[2]);
        const Matrix Head2IPMat(argv[4]);
        const Matrix SourceMat(argv[3]);

        const Matrix& InternalPotGainMat = (Head2IPMat*HeadMatInv)*SourceMat;

        InternalPotGainMat.save(argv[5]);

    } else {

        error(argv[0],true);
    }

    // Stop Chrono

    const auto end_time = std::chrono::system_clock::now();
    dispEllapsed(end_time-start_time);

    return 0;
}

void
getHelp(char** argv) {

    std::cout << argv[0] <<" [-option] [filepaths...]" << std::endl << std::endl;

    std::cout << "-option :" << std::endl;
    std::cout << "   -EEG :   Compute the gain for EEG " << std::endl;
    std::cout << "            Filepaths are in order :" << std::endl;
    std::cout << "            HeadMatInv, SourceMat, Head2EEGMat, EEGGainMatrix" << std::endl;
    std::cout << "            bin Matrix" << std::endl << std::endl;

    std::cout << "   -MEG :   Compute the gain for MEG " << std::endl;
    std::cout << "            Filepaths are in order :" << std::endl;
    std::cout << "            HeadMatInv, SourceMat, Head2MEGMat, Source2MEGMat, MEGGainMatrix" << std::endl;
    std::cout << "            Matrix (.bin or .txt)" << std::endl << std::endl;

    std::cout << "   -InternalPotential or -IP :   Compute the gain for internal potentials, measured within the volume " << std::endl;
    std::cout << "            Filepaths are in order :" << std::endl;
    std::cout << "            HeadMatInv, SourceMat, Head2IPMat, Source2IPMat" << std::endl;
    std::cout << "            InternalPotential gain Matrix (.txt)" << std::endl << std::endl;

    std::cout << "   -EITInternalPotential or -EITIP :   Compute the gain for internal potentials using boundary normal current as input" << std::endl;
    std::cout << "            Filepaths are in order :" << std::endl;
    std::cout << "            HeadMatInv, SourceMat, Head2IPMat" << std::endl;
    std::cout << "            EITInternalPotential gain Matrix" << std::endl << std::endl;

    std::cout << "   -EEGadjoint :   Compute the gain for EEG " << std::endl;
    std::cout << "            Filepaths are in order :" << std::endl;
    std::cout << "            geometry file (.geom)" << std::endl;
    std::cout << "            conductivity file (.cond)" << std::endl;
    std::cout << "            dipoles positions and orientations" << std::endl;
    std::cout << "            HeadMat, Head2EEGMat, EEGGainMatrix" << std::endl;
    std::cout << "            bin Matrix" << std::endl << std::endl;

    std::cout << "   -MEGadjoint :   Compute the gain for MEG " << std::endl;
    std::cout << "            Filepaths are in order :" << std::endl;
    std::cout << "            geometry file (.geom)" << std::endl;
    std::cout << "            conductivity file (.cond)" << std::endl;
    std::cout << "            dipoles positions and orientations" << std::endl;
    std::cout << "            HeadMat, Head2MEGMat, Source2MEGMat, MEGGainMatrix" << std::endl;
    std::cout << "            bin Matrix" << std::endl << std::endl;

    std::cout << "   -EEGMEGadjoint :   Compute the gain for EEG and MEG at the same time" << std::endl;
    std::cout << "            Filepaths are in order :" << std::endl;
    std::cout << "            geometry file (.geom)" << std::endl;
    std::cout << "            conductivity file (.cond)" << std::endl;
    std::cout << "            dipoles positions and orientations" << std::endl;
    std::cout << "            HeadMat, Head2EEGMat, Head2MEGMat, Source2MEGMat, EEGGainMatrix, MEGGainMatrix" << std::endl;
    std::cout << "            bin Matrix" << std::endl << std::endl;

    exit(0);
}
