/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Maureen.Clerc.AT.sophia.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.sophia.inria.fr)

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

#include "matrix.h"
#include "matrix.h"
#include "symmatrix.h"
#include "vector.h"
#include "cpuChrono.h"
#include "gain.h"

using namespace std;
using namespace OpenMEEG;

void getHelp(char** argv);

int main(int argc, char **argv)
{
    print_version(argv[0]);

    if(argc<2)
    {
        cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
        return 0;
    }

    if ((!strcmp(argv[1],"-h")) | (!strcmp(argv[1],"--help"))) getHelp(argv);

    // Start Chrono
    cpuChrono C;
    C.start();

    disp_argv(argc,argv);

    // declaration of argument variables
    string Option=string(argv[1]);
    if(argc<5)
    {
        cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
        return 0;
    }

    // for use with EEG DATA
    if(!strcmp(argv[1],"-EEG"))
    {
        if(argc<6)
        {
            cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
            return 0;
        }

        Matrix EEGGainMatrix;

        { // Avoiding to store all matrices at the same time
            LinOpInfo matinfo = OpenMEEG::maths::info(argv[3]);
            SymMatrix HeadMatInv;
            HeadMatInv.load(argv[2]);
            EEGGainMatrix = HeadMatInv(0,HeadMatInv.nlin()-1,0,matinfo.nlin()-1); // reducedHeadInvMatrix
        }
        {
            SparseMatrix Head2EEGMatrix;
            Head2EEGMatrix.load(argv[4]);
            EEGGainMatrix = Head2EEGMatrix*EEGGainMatrix;
        }
        {
            Matrix SourceMatrix;
            SourceMatrix.load(argv[3]);
            EEGGainMatrix = EEGGainMatrix*SourceMatrix;
        }

        EEGGainMatrix.save(argv[5]);
    }
    // for use with MEG DATA
    else if(!strcmp(argv[1],"-MEG"))
    {
        if(argc<7)
        {
            cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
            return 0;
        }
        Matrix MEGGainMatrix;

        { // Avoiding to store all matrices at the same time
            LinOpInfo matinfo = OpenMEEG::maths::info(argv[3]);
            SymMatrix HeadMatInv;
            HeadMatInv.load(argv[2]);
            MEGGainMatrix = HeadMatInv(0,HeadMatInv.nlin()-1,0,matinfo.nlin()-1); // reducedLhsInvMatrix
        }
        {
            Matrix Head2MEGMatrix;
            Head2MEGMatrix.load(argv[4]);
            MEGGainMatrix = Head2MEGMatrix*MEGGainMatrix;
        }
        {
            Matrix SourceMatrix;
            SourceMatrix.load(argv[3]);
            MEGGainMatrix = MEGGainMatrix*SourceMatrix;
        }
        {
            Matrix Source2MEGMatrix;
            Source2MEGMatrix.load(argv[5]);
            MEGGainMatrix += Source2MEGMatrix;
        }

        MEGGainMatrix.save(argv[6]);
    }
    else if(!strcmp(argv[1],"-VolEEG"))
    {
        if(argc<6)
        {
            cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
            return 0;
        }
        Matrix VolEEGGainMatrix;

        { // Avoiding to store all matrices at the same time
            Matrix Surf2Vol;
            Surf2Vol.load(argv[4]); 
            SymMatrix HeadMatInv;
            HeadMatInv.load(argv[2]);
            VolEEGGainMatrix = Surf2Vol*HeadMatInv(0,Surf2Vol.ncol()-1,0,HeadMatInv.ncol()-1);
        }
        {
            Matrix EITStim;
            EITStim.load(argv[3]);
            VolEEGGainMatrix = VolEEGGainMatrix*EITStim;
        }
        VolEEGGainMatrix.save(argv[5]);
    }
    else
    {
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
    cout << "            HeadMatInv, SourceMat, Head2MEGMatrix, Source2MEGMatrix, MEGGainMatrix" << endl;
    cout << "            Matrix (.bin or .txt)" << endl << endl;

    cout << "   -VolEEG :   Compute the gain for EEG, measured within the volume " << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            inputs: HeadMatInv, SourceMat, Surf2VolMat" << endl;
    cout << "            output: VolEEGgain Matrix (.txt)" << endl << endl;

    exit(0);
}
