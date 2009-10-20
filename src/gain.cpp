/* FILE: $Id$ */

/*
Project Name : OpenMEEG

version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

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
    om_print_version(argv);

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
            size_t nn,mm;
            Matrix::readDimsBin(argv[3],nn,mm); // read nb lines of SourceMatrix without allocating it
            SymMatrix HeadMatInv;
            HeadMatInv.loadBin(argv[2]);
            EEGGainMatrix = HeadMatInv(0,HeadMatInv.nlin()-1,0,nn-1); // reducedHeadInvMatrix
        }
        {
            SparseMatrix Head2EEGMatrix;
            Head2EEGMatrix.loadBin(argv[4]);
            EEGGainMatrix = Head2EEGMatrix*EEGGainMatrix;
        }
        {
            Matrix SourceMatrix;
            SourceMatrix.loadBin(argv[3]);
            EEGGainMatrix = EEGGainMatrix*SourceMatrix;
        }

        EEGGainMatrix.saveBin(argv[5]);
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
            size_t nn,mm;
            Matrix::readDimsBin(argv[3],nn,mm); // read nb lines of SourceMatrix without allocating it
            SymMatrix HeadMatInv;
            HeadMatInv.loadBin(argv[2]);
            MEGGainMatrix = HeadMatInv(0,HeadMatInv.nlin()-1,0,nn-1); // reducedLhsInvMatrix
        }
        {
            Matrix Head2MEGMatrix;
            Head2MEGMatrix.loadBin(argv[4]);
            MEGGainMatrix = Head2MEGMatrix*MEGGainMatrix;
        }
        {
            Matrix SourceMatrix;
            SourceMatrix.loadBin(argv[3]);
            MEGGainMatrix = MEGGainMatrix*SourceMatrix;
        }
        {
            Matrix Source2MEGMatrix;
            Source2MEGMatrix.loadBin(argv[5]);
            MEGGainMatrix += Source2MEGMatrix;
        }

        MEGGainMatrix.saveBin(argv[6]);
    }
    else if(!strcmp(argv[1],"-VolPotEIT"))
    {
        if(argc<6)
        {
            cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
            return 0;
        }
        Matrix VolPotEITGainMatrix;

        { // Avoiding to store all matrices at the same time
            Matrix Surf2Vol;
            Surf2Vol.loadBin(argv[2]); 
            SymMatrix HeadMatInv;
            HeadMatInv.loadBin(argv[3]);
            VolPotEITGainMatrix = Surf2Vol*HeadMatInv(0,Surf2Vol.ncol()-1,0,HeadMatInv.ncol()-1);
        }
        {
            Matrix EITStim;
            EITStim.loadBin(argv[4]);
            VolPotEITGainMatrix = VolPotEITGainMatrix*EITStim;
        }
        VolPotEITGainMatrix.saveTxt(argv[5]);
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

    cout << "   -VolPotEIT :   Compute the gain for EIT, measured within the volume " << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            inputs: Surf2VolMat, HeadMatInv, EITStimMatrix," << endl;
    cout << "            output: VolPotEITgain Matrix (.txt)" << endl << endl;

    exit(0);

}
