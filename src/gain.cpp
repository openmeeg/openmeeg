#include <cstring>

#include "matrice.h"
#include "symmatrice.h"
#include "vecteur.h"
#include "cpuChrono.h"

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

    cout << endl << "| ------ " << argv[0] << " -------" << endl;
    for( int i = 1; i < argc; i += 1 )
    {
        cout << "| " << argv[i] << endl;
    }
    cout << "| -----------------------" << endl;

    // declaration of argument variables
    string Option=string(argv[1]);
    if(argc<5)
    {
        cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
        return 0; 
    }
    symmatrice LhsInvMatrix; 
    matrice RhsMatrix; 
    matrice V2MegMatrix;
    matrice S2MegMatrix;
    matrice MegGainMatrix;
    matrice V2EegMatrix;
    matrice EegGainMatrix;

    matrice reducedLhsInvMatrix; 

    // for use with EEG DATA
    if(!strcmp(argv[1],"-EEG"))
    {
        if(argc<6)
        {
            cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
            return 0; 
        }
        LhsInvMatrix.loadBin(argv[2]);
        RhsMatrix.loadBin(argv[3]);
        reducedLhsInvMatrix=matrice(LhsInvMatrix)(0,LhsInvMatrix.nlin()-1,0,RhsMatrix.nlin()-1);

        V2EegMatrix.loadBin(argv[4]);
        EegGainMatrix=(V2EegMatrix*reducedLhsInvMatrix)*RhsMatrix;
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
        LhsInvMatrix.loadBin(argv[2]);
        RhsMatrix.loadBin(argv[3]);
        reducedLhsInvMatrix=matrice(LhsInvMatrix)(0,LhsInvMatrix.nlin()-1,0,RhsMatrix.nlin()-1);

        V2MegMatrix.loadBin(argv[4]);
        S2MegMatrix.loadBin(argv[5]);
        MegGainMatrix=S2MegMatrix+(V2MegMatrix*reducedLhsInvMatrix)*RhsMatrix;
        MegGainMatrix.saveBin(argv[6]);
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

    exit(0);

}


