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

#include <cstring>

#include "matrice.h"
#include "symmatrice.h"
#include "vecteur.h"
#include "cpuChrono.h"

using namespace std;

void getHelp(char** argv);

int main(int argc, char **argv)
{
    if(argc==1)
    {
        cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
        return 0;
    }

    if ((!strcmp(argv[1],"-h")) | (!strcmp(argv[1],"--help"))) getHelp(argv);

    disp_argv(argc,argv);

    // Start Chrono
    cpuChrono C;
    C.start();

    symmatrice LhsMatrix;
    symmatrice LhsInvMatrix;

    LhsMatrix.loadBin(argv[1]);
    LhsInvMatrix=LhsMatrix.inverse();
    LhsInvMatrix.saveBin(argv[2]);

    // Stop Chrono
    C.stop();
    C.dispEllapsed();

    return 0;
}

void getHelp(char** argv)
{
    cout << argv[0] <<" [-option] [filepaths...]" << endl << endl;

    cout << "   Inverse LHS " << endl;
    cout << "   Filepaths are in order :" << endl;
    cout << "       LhsMatrix (bin), LhsInvMatrix (bin)" << endl << endl;

    exit(0);
}

