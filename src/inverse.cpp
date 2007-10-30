#include "cpuChrono.h"
#include "om_utils.h"
#include "inverse.h"

using namespace std;

void getHelp(char** argv);

int main(int argc, char **argv)
{
    // Makefile for global usage
    if(argc==1)
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
    matrice GainMatrix;
    sparse_matrice SmoothMatrix;
    vecteur AiVector;
    matrice Data;
    matrice EstimatedSourcesData;
    double SmoothWeight;
    string SmoothType;

    GainMatrix.loadBin(argv[1]);
    SmoothMatrix.loadBin(argv[2]);
    AiVector.loadBin(argv[3]);
    Data.loadTxt(argv[4]);
    SmoothWeight = atof(argv[6] );
    SmoothType   = string(argv[7]);

    bool Heat = SmoothType==string("HEAT");
    bool Mn   = SmoothType==string("MN");
    bool Tv   = SmoothType==string("TV");

    if (!Tv && !Mn && !Heat) {
        std::cerr << "Unknown Smoothtype :  " << SmoothType << std::endl;
        std::cerr << "Should be HEAT , MN or TV" << std::endl;
        exit(1);
    }

    if(Tv)
    {
        size_t MaxNbIter   = (size_t) atoi(argv[8]);
        double StoppingTol = atof(argv[9]);
    
        TV_inverse_matrice EstimatedSourcesData(Data,GainMatrix,SmoothMatrix,AiVector,SmoothWeight,MaxNbIter,StoppingTol);
        EstimatedSourcesData.saveTxt(argv[5]);
    }

    if(Mn)
    {
        MN_inverse_matrice EstimatedSourcesData(Data,GainMatrix,SmoothWeight);
        EstimatedSourcesData.saveTxt(argv[5]);
    }

    if(Heat)
    {
        HEAT_inverse_matrice EstimatedSourcesData(Data,GainMatrix,SmoothMatrix,SmoothWeight);
        EstimatedSourcesData.saveTxt(argv[5]);
    }

    // Stop Chrono
    C.stop();
    C.dispEllapsed();

    return 0;
}

void getHelp(char** argv)
{
    cout << argv[0] <<"[filepaths...]" << endl << endl;
    cout << "Compute the inverse for MEG/EEG " << endl;
    cout << "\tFilepaths are in order :" << endl;
    cout << "\tGainMatrix, SmoothMatrix, AiVector, RealData, EstimatedSourcesData, SmoothWeight, SmoothType, MaxNbIter, StoppingTol" << endl << endl;
    exit(0);
}
