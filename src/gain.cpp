/* FILE: $Id$ */

/*
Project Name : OpenMEEG

author            : $Author$
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

#include "matrice.h"
#include "matrice_dcl.h"
#include "symmatrice.h"
#include "vecteur.h"
#include "cpuChrono.h"
#include "gain.h"

using namespace std;

void getHelp(char** argv);

int main(int argc, char **argv)
{

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

        matrice EegGainMatrix;

        { // Avoiding to store all matrices at the same time
            size_t nn,mm;
            matrice::readDimsBin(argv[3],nn,mm); // read nb lines of RhsMatrix without allocating it
            symmatrice LhsInvMatrix;
            LhsInvMatrix.loadBin(argv[2]);
            EegGainMatrix = matrice(LhsInvMatrix)(0,LhsInvMatrix.nlin()-1,0,nn-1); // reducedLhsInvMatrix
        }
        {
            matrice V2EegMatrix;
            V2EegMatrix.loadBin(argv[4]);
            EegGainMatrix = V2EegMatrix*EegGainMatrix;
        }
        {
            matrice RhsMatrix;
            RhsMatrix.loadBin(argv[3]);
            EegGainMatrix = EegGainMatrix*RhsMatrix;
        }

        EegGainMatrix.saveBin(argv[5]);
    }
    // for use with MEG DATA
    else if(!strcmp(argv[1],"-MEG"))
    {
        if(argc<7)
        {
            cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
            return 0;
        }
        matrice MegGainMatrix;

        { // Avoiding to store all matrices at the same time
            size_t nn,mm;
            matrice::readDimsBin(argv[3],nn,mm); // read nb lines of RhsMatrix without allocating it
            symmatrice LhsInvMatrix;
            LhsInvMatrix.loadBin(argv[2]);
            MegGainMatrix = matrice(LhsInvMatrix)(0,LhsInvMatrix.nlin()-1,0,nn-1); // reducedLhsInvMatrix
        }
        {
            matrice V2MegMatrix;
            V2MegMatrix.loadBin(argv[4]);
            MegGainMatrix = V2MegMatrix*MegGainMatrix;
        }
        {
            matrice RhsMatrix;
            RhsMatrix.loadBin(argv[3]);
            MegGainMatrix = MegGainMatrix*RhsMatrix;
        }
        {
            matrice S2MegMatrix;
            S2MegMatrix.loadBin(argv[5]);
            MegGainMatrix += S2MegMatrix;
        }

        MegGainMatrix.saveBin(argv[6]);
    }
    else if(!strcmp(argv[1],"-VolPotEIT"))
    {
        if(argc<6)
        {
            cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
            return 0;
        }
        matrice VolPotEITGainMatrix;

        { // Avoiding to store all matrices at the same time
	  matrice SurfToVol;
	  SurfToVol.loadBin(argv[2]); 
	  symmatrice LhsInvMatrix;
	  LhsInvMatrix.loadBin(argv[3]);
	  VolPotEITGainMatrix = SurfToVol*matrice(LhsInvMatrix)(0,SurfToVol.ncol()-1,0,LhsInvMatrix.ncol()-1);
        }
        {
	  matrice EITStim;
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
    cout << "            LhsInvMatrix, RhsMatrix, V2EegMatrix, EegGainMatrix" << endl;
    cout << "            bin matrix" << endl << endl;

    cout << "   -MEG :   Compute the gain for MEG " << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            LhsInvMatrix, RhsMatrix, V2MegMatrix, S2MegMatrix, MegGainMatrix" << endl;
    cout << "            matrix (.bin or .txt)" << endl << endl;

    cout << "   -VolPotEIT :   Compute the gain for EIT, measured within the volume " << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            inputs: SurfToVolMatrix, LhsInvMatrix, EITStimMatrix," << endl;
    cout << "            output: VolPotEITgain matrix (.txt)" << endl << endl;

    exit(0);

}


