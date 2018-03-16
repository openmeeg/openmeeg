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

#include <cstring>
#include <cpuChrono.h>
#include <gain.h>

using namespace std;
using namespace OpenMEEG;

void getHelp(char** argv);

int main(int argc, char **argv)
{
    print_version(argv[0]);

    if ( argc<2 ) {
        cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
        return 0;
    }

    if ( (!strcmp(argv[1], "-h")) | (!strcmp(argv[1], "--help")) ) {
        getHelp(argv);
    }

    // Start Chrono
    cpuChrono C;
    C.start();

    disp_argv(argc, argv);

    // declaration of argument variables
    string Option=string(argv[1]);
    if ( argc<5 ) {
        cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
        return 0;
    }

    // for use with EEG DATA
    if ( !strcmp(argv[1], "-EEG") ) {
        if ( argc<6 ) {
            cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
            return 0;
        }
        //  We split the 2 matrix multiplications in order to spare memory.
        //  This is also why we do not use GainEEG...

        const SymMatrix HeadMatInv(argv[2]);
        const SparseMatrix Head2EEGMat(argv[4]);
        const Matrix& tmp = Head2EEGMat*HeadMatInv;
        const Matrix SourceMat(argv[3]);
        const Matrix& EEGGainMat = tmp*SourceMat;
        EEGGainMat.save(argv[5]);
    }
    // compute the gain matrix with the adjoint method for use with EEG DATA
    else if ( !strcmp(argv[1], "-EEGadjoint") ) {
        if ( argc<8 ) {
            cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
            return 0;
        }
        Geometry geo;
        geo.read(argv[2], argv[3]);
        const Matrix dipoles(argv[4]);
        const SymMatrix HeadMat(argv[5]);
        const SparseMatrix Head2EEGMat(argv[6]);

        const GainEEGadjoint EEGGainMat(geo, dipoles, HeadMat, Head2EEGMat);
        EEGGainMat.save(argv[7]);
    }
    // for use with MEG DATA
    else if ( !strcmp(argv[1], "-MEG") ) {
        if ( argc<7 ) {
            cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
            return 0;
        }
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
    }
    // compute the gain matrix with the adjoint method for use with MEG DATA
    else if ( !strcmp(argv[1], "-MEGadjoint") ) {
        if ( argc<9 ) {
            cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
            return 0;
        }
        Geometry geo;
        geo.read(argv[2], argv[3]);
        const Matrix dipoles(argv[4]);
        const SymMatrix HeadMat(argv[5]);
        const Matrix Head2MEGMat(argv[6]);
        const Matrix Source2MEGMat(argv[7]);

        const GainMEGadjoint MEGGainMat(geo, dipoles, HeadMat, Head2MEGMat, Source2MEGMat);
        MEGGainMat.save(argv[8]);
    }
    // compute the gain matrices with the adjoint method for use with EEG and MEG DATA
    else if ( !strcmp(argv[1], "-EEGMEGadjoint") ) {
        if ( argc<11 ) {
            cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
            return 0;
        }
        Geometry geo;
        geo.read(argv[2], argv[3]);
        const Matrix dipoles(argv[4]);
        const SymMatrix HeadMat(argv[5]);
        const SparseMatrix Head2EEGMat(argv[6]);
        const Matrix Head2MEGMat(argv[7]);
        const Matrix Source2MEGMat(argv[8]);

        const GainEEGMEGadjoint EEGMEGGainMat(geo, dipoles, HeadMat, Head2EEGMat, Head2MEGMat, Source2MEGMat);
        EEGMEGGainMat.saveEEG(argv[9]);
        EEGMEGGainMat.saveMEG(argv[10]);
    }
    else if ( (!strcmp(argv[1], "-InternalPotential"))|(!strcmp(argv[1], "-IP")) ) {
        if ( argc<7 ) {
            cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
            return 0;
        }
        const SymMatrix HeadMatInv(argv[2]);
        const Matrix Head2IPMat(argv[4]);

        const Matrix& tmp1 = Head2IPMat*HeadMatInv;
        const Matrix SourceMat(argv[3]);
        const Matrix& tmp2 = tmp1*SourceMat;
        const Matrix Source2IPMat(argv[5]);

        const Matrix& InternalPotGainMat = Source2IPMat+tmp2;
        InternalPotGainMat.save(argv[6]);
    }
    else if ( (!strcmp(argv[1], "-EITInternalPotential"))||(!strcmp(argv[1], "-EITIP")) ) {
        if ( argc<6 ){
            cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
            return 0;
        }
        const SymMatrix HeadMatInv(argv[2]);
        const Matrix Head2IPMat(argv[4]);

        const Matrix& tmp = Head2IPMat*HeadMatInv;

        const Matrix SourceMat(argv[3]);

        const Matrix& InternalPotGainMat = tmp*SourceMat;
        InternalPotGainMat.save(argv[5]);
    }
    else {
        cerr << "Error: unknown option. \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
    }

    // Stop Chrono
    C.stop();
    C.dispEllapsed();

    return 0;
}

void getHelp(char** argv)
{
    cout << argv[0] <<" [-option] [filepaths...]" << endl << endl;

    cout << "-option :" << endl;
    cout << "   -EEG :   Compute the gain for EEG " << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            HeadMatInv, SourceMat, Head2EEGMat, EEGGainMatrix" << endl;
    cout << "            bin Matrix" << endl << endl;

    cout << "   -MEG :   Compute the gain for MEG " << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            HeadMatInv, SourceMat, Head2MEGMat, Source2MEGMat, MEGGainMatrix" << endl;
    cout << "            Matrix (.bin or .txt)" << endl << endl;

    cout << "   -InternalPotential or -IP :   Compute the gain for internal potentials, measured within the volume " << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            HeadMatInv, SourceMat, Head2IPMat, Source2IPMat" << endl;
    cout << "            InternalPotential gain Matrix (.txt)" << endl << endl;

    cout << "   -EITInternalPotential or -EITIP :   Compute the gain for internal potentials using boundary normal current as input" << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            HeadMatInv, SourceMat, Head2IPMat" << endl;
    cout << "            EITInternalPotential gain Matrix" << endl << endl;

    cout << "   -EEGadjoint :   Compute the gain for EEG " << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            geometry file (.geom)" << endl;
    cout << "            conductivity file (.cond)" << endl;
    cout << "            dipoles positions and orientations" << endl;
    cout << "            HeadMat, Head2EEGMat, EEGGainMatrix" << endl;
    cout << "            bin Matrix" << endl << endl;

    cout << "   -MEGadjoint :   Compute the gain for MEG " << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            geometry file (.geom)" << endl;
    cout << "            conductivity file (.cond)" << endl;
    cout << "            dipoles positions and orientations" << endl;
    cout << "            HeadMat, Head2MEGMat, Source2MEGMat, MEGGainMatrix" << endl;
    cout << "            bin Matrix" << endl << endl;

    cout << "   -EEGMEGadjoint :   Compute the gain for EEG and MEG at the same time" << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            geometry file (.geom)" << endl;
    cout << "            conductivity file (.cond)" << endl;
    cout << "            dipoles positions and orientations" << endl;
    cout << "            HeadMat, Head2EEGMat, Head2MEGMat, Source2MEGMat, EEGGainMatrix, MEGGainMatrix" << endl;
    cout << "            bin Matrix" << endl << endl;

    exit(0);
}
