/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
lastrevision      : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
*/

#include "cpuChrono.h"
#include "om_utils.h"
#include "forward.h"

using namespace std;

void getHelp(char** argv);

int main(int argc, char **argv)
{
    if(argc==1)
    {
        cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
        return 0;
    }

    if ((!strcmp(argv[1],"-h")) | (!strcmp(argv[1],"--help"))) {
        getHelp(argv);
        return 0;
    }

    if(argc < 5)
    {
        cerr << "Bad arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
        exit(1);
    }

    // Start Chrono
    cpuChrono C;
    C.start();

    disp_argv(argc,argv);

    // declaration of argument variables======================================================================
    matrice GainMatrix;
    matrice RealSourcesData;
    double NoiseLevel;

    GainMatrix.loadBin(argv[1]);
    RealSourcesData.loadTxt(argv[2]);
    NoiseLevel = atof(argv[4]);

    Forward_matrice SimulatedData(GainMatrix,RealSourcesData,NoiseLevel);

    // write output variables ===================================================================================
    SimulatedData.saveTxt(argv[3]);

    // Stop Chrono
    C.stop();
    C.dispEllapsed();

    return 0;
}

void getHelp(char** argv)
{
    cout << argv[0] << " [filepaths...]" << endl << endl;

    cout << "   Compute the forward problem " << endl;
    cout << "   Filepaths are in order :" << endl;
    cout << "   GainMatrix (bin), RealSourcesData (txt), SimulatedData (txt), NoiseLevel (float)" << endl << endl;

    exit(0);
}


