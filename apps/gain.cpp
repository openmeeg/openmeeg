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

void help(const char* cmd_name);

inline void
error(const char* command,const bool unknown_option=false) {
    std::cerr << "Error: " << ((unknown_option) ? "Unknown option." : "Not enough arguments.") << std::endl
              << "Please try \"" << command << " -h\" or \"" << command << " --help \" \n" << std::endl;
    exit(1);
}

int
main(int argc,char** argv) {

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

    constexpr char geomfileopt[]       = "geometry file";
    constexpr char condfileopt[]       = "conductivity file";
    constexpr char outputfileopt[]     = "output file";
    constexpr char dipolefileopt[]     = "dipoles file";
    constexpr char HMfileopt[]         = "head matrix";
    constexpr char invHMfileopt[]      = "inverse head matrix";
    constexpr char sourcematfileopt[]  = "source matrix";
    constexpr char h2eegfileopt[]      = "head to EEG matrix";
    constexpr char h2megfileopt[]      = "head to MEG matrix";
    constexpr char h2ipfileopt[]       = "head to internal potential matrix";
    constexpr char source2megfileopt[] = "sources to MEG matrix";
    constexpr char source2ipfileopt[]  = "sources to internal potential matrix";

    // Start Chrono

    const auto start_time = std::chrono::system_clock::now();

    unsigned num_options = 0;

    const auto& EEGparms = {
        invHMfileopt, sourcematfileopt, h2eegfileopt, outputfileopt
    };

    if (char** opt_parms = cmd.option("-EEG",EEGparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        // EEG DATA

        //  Split the 2 matrix multiplications in order to spare memory.
        //  This is why we do not use GainEEG...

        const SymMatrix    HeadMatInv(opt_parms[1]);
        const SparseMatrix Head2EEGMat(opt_parms[3]);
        const Matrix& tmp = Head2EEGMat*HeadMatInv;
        const Matrix SourceMat(opt_parms[2]);
        const Matrix& EEGGainMat = tmp*SourceMat;
        EEGGainMat.save(opt_parms[4]);
    }

    const auto& EEGAdjointparms = {
        geomfileopt, condfileopt, dipolefileopt, HMfileopt, h2eegfileopt, outputfileopt
    };

    if (char** opt_parms = cmd.option("-EEGadjoint",EEGAdjointparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        // Compute gain matrix with the adjoint method for use with EEG DATA

        Geometry geo(opt_parms[1],opt_parms[2]);
        const Matrix dipoles(opt_parms[3]);
        const SymMatrix HeadMat(opt_parms[4]);
        const SparseMatrix Head2EEGMat(opt_parms[5]);

        const GainEEGadjoint EEGGainMat(geo, dipoles, HeadMat, Head2EEGMat);
        EEGGainMat.save(opt_parms[6]);
    }

    const auto& MEGparms = {
        invHMfileopt, sourcematfileopt, h2megfileopt, source2megfileopt, outputfileopt
    };

    if (char** opt_parms = cmd.option("-MEG",MEGparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        // MEG DATA

        //  We split the 3 matrix multiplications in order to spare memory.
        //  This is also why we do not use GainMEG...

        const SymMatrix HeadMatInv(opt_parms[1]);
        const Matrix Head2MEGMat(opt_parms[3]);
        const Matrix& tmp1 = Head2MEGMat*HeadMatInv;
        const Matrix SourceMat(opt_parms[2]);
        const Matrix& tmp2 = tmp1*SourceMat;
        const Matrix Source2MEGMat(opt_parms[4]);
        const Matrix MEGGainMat = Source2MEGMat+tmp2;
        MEGGainMat.save(opt_parms[5]);
    }

    const auto& MEGAdjointparms = {
        geomfileopt, condfileopt, dipolefileopt, HMfileopt, h2megfileopt, source2megfileopt, outputfileopt
    };

    if (char** opt_parms = cmd.option("-MEGadjoint",MEGAdjointparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        // Compute the gain matrix with the adjoint method for use with MEG DATA

        Geometry geo(opt_parms[1],opt_parms[2]);
        const Matrix dipoles(opt_parms[3]);
        const SymMatrix HeadMat(opt_parms[4]);
        const Matrix Head2MEGMat(opt_parms[5]);
        const Matrix Source2MEGMat(opt_parms[6]);

        const GainMEGadjoint MEGGainMat(geo, dipoles, HeadMat, Head2MEGMat, Source2MEGMat);
        MEGGainMat.save(opt_parms[7]);
    }

    const auto& EEGMEGAdjointparms = {
        geomfileopt, condfileopt, dipolefileopt, HMfileopt, h2eegfileopt, h2megfileopt, source2megfileopt,
        "EEG output file", "MEG output file"
    };

    if (char** opt_parms = cmd.option("-EEGMEGadjoint",EEGMEGAdjointparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        // Compute the gain matrices with the adjoint method for EEG and MEG DATA

        Geometry geo(opt_parms[1],opt_parms[2]);
        const Matrix dipoles(opt_parms[3]);
        const SymMatrix HeadMat(opt_parms[4]);
        const SparseMatrix Head2EEGMat(opt_parms[5]);
        const Matrix Head2MEGMat(opt_parms[6]);
        const Matrix Source2MEGMat(opt_parms[7]);

        const GainEEGMEGadjoint EEGMEGGainMat(geo,dipoles,HeadMat,Head2EEGMat,Head2MEGMat,Source2MEGMat);
        EEGMEGGainMat.saveEEG(opt_parms[8]);
        EEGMEGGainMat.saveMEG(opt_parms[9]);
    }

    const auto& IPparms = {
        invHMfileopt, sourcematfileopt, h2ipfileopt, source2ipfileopt, outputfileopt
    };

    if (char** opt_parms = cmd.option({ "-InternalPotential", "-IP", "-ip" },IPparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        const SymMatrix HeadMatInv(opt_parms[1]);
        const Matrix Head2IPMat(opt_parms[3]);

        const Matrix& tmp1 = Head2IPMat*HeadMatInv;
        const Matrix SourceMat(opt_parms[2]);
        const Matrix& tmp2 = tmp1*SourceMat;
        const Matrix Source2IPMat(opt_parms[4]);

        const Matrix& InternalPotGainMat = Source2IPMat+tmp2;
        InternalPotGainMat.save(opt_parms[5]);
    }

    const auto& EITIPparms = {
        invHMfileopt, sourcematfileopt, h2ipfileopt, outputfileopt
    };

    if (char** opt_parms = cmd.option({ "-EITInternalPotential", "-EITIP", "-eitip" },EITIPparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        const SymMatrix HeadMatInv(opt_parms[1]);
        const Matrix SourceMat(opt_parms[2]);
        const Matrix Head2IPMat(opt_parms[3]);

        const Matrix& InternalPotGainMat = (Head2IPMat*HeadMatInv)*SourceMat;

        InternalPotGainMat.save(opt_parms[4]);
    }

    if (num_options==0) {
        std::cerr << "Unknown argument: " << argv[1] << std::endl;
        exit(1);
    }

    // Stop Chrono

    const auto end_time = std::chrono::system_clock::now();
    dispEllapsed(end_time-start_time);

    return 0;
}

void
help(const char* cmd_name) {

    std::cout << cmd_name << " [-option] [filepaths...]" << std::endl << std::endl;

    std::cout << "-option :" << std::endl
              << "   -EEG :   Compute the gain for EEG " << std::endl
              << "            Filepaths are in order :" << std::endl
              << "            HeadMatInv, SourceMat, Head2EEGMat, EEGGainMatrix" << std::endl
              << "            bin Matrix" << std::endl << std::endl;

    std::cout << "   -MEG :   Compute the gain for MEG " << std::endl
              << "            Filepaths are in order :" << std::endl
              << "            HeadMatInv, SourceMat, Head2MEGMat, Source2MEGMat, MEGGainMatrix" << std::endl
              << "            Matrix (.bin or .txt)" << std::endl << std::endl;

    std::cout << "   -InternalPotential or -IP :   Compute the gain for internal potentials, measured within the volume " << std::endl
              << "            Filepaths are in order :" << std::endl
              << "            HeadMatInv, SourceMat, Head2IPMat, Source2IPMat" << std::endl
              << "            InternalPotential gain Matrix (.txt)" << std::endl << std::endl;

    std::cout << "   -EITInternalPotential or -EITIP :   Compute the gain for internal potentials using boundary normal current as input" << std::endl
              << "            Filepaths are in order :" << std::endl
              << "            HeadMatInv, SourceMat, Head2IPMat" << std::endl
              << "            EITInternalPotential gain Matrix" << std::endl << std::endl;

    std::cout << "   -EEGadjoint :   Compute the gain for EEG " << std::endl
              << "            Filepaths are in order :" << std::endl
              << "            geometry file (.geom)" << std::endl
              << "            conductivity file (.cond)" << std::endl
              << "            dipoles positions and orientations" << std::endl
              << "            HeadMat, Head2EEGMat, EEGGainMatrix" << std::endl
              << "            bin Matrix" << std::endl << std::endl;

    std::cout << "   -MEGadjoint :   Compute the gain for MEG " << std::endl
              << "            Filepaths are in order :" << std::endl
              << "            geometry file (.geom)" << std::endl
              << "            conductivity file (.cond)" << std::endl
              << "            dipoles positions and orientations" << std::endl
              << "            HeadMat, Head2MEGMat, Source2MEGMat, MEGGainMatrix" << std::endl
              << "            bin Matrix" << std::endl << std::endl;

    std::cout << "   -EEGMEGadjoint :   Compute the gain for EEG and MEG at the same time" << std::endl
              << "            Filepaths are in order :" << std::endl
              << "            geometry file (.geom)" << std::endl
              << "            conductivity file (.cond)" << std::endl
              << "            dipoles positions and orientations" << std::endl
              << "            HeadMat, Head2EEGMat, Head2MEGMat, Source2MEGMat, EEGGainMatrix, MEGGainMatrix" << std::endl
              << "            bin Matrix" << std::endl << std::endl;
}
