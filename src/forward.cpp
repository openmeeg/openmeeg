#include "matrice.h"
#include "symmatrice.h"
#include "vecteur.h"
#include "cpuChrono.h"
#include "om_utils.h"

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

    // Start Chrono
    cpuChrono C;
    C.start();

    cout << endl << "| ------ " << argv[0] << " -------" << endl;
    for( int i = 1; i < argc; i += 1 )
    {
        cout << "| " << argv[i] << endl;
    }
    cout << "| -----------------------" << endl;

    // declaration of argument variables======================================================================
    string Option;
    matrice MegGainMatrix;
    matrice EegGainMatrix;
    matrice RealSourcesData;
    matrice SimulatedMegData;
    double MegNoiseLevel;
    matrice SimulatedEegData;
    double EegNoiseLevel;
    
    // for use with EEG DATA
    if(!strcmp(argv[1],"-EEG"))
    {
        Option=string(argv[1]);
        EegGainMatrix.loadBin(argv[2]);
        RealSourcesData.loadTxt(argv[3]);
        EegNoiseLevel=atof(argv[5]);
    }
    // for use with MEG DATA
    else if(!strcmp(argv[1],"-MEG"))
    {
        Option=string(argv[1]);
        MegGainMatrix.loadBin(argv[2]);
        RealSourcesData.loadTxt(argv[3]);
        MegNoiseLevel=atof(argv[5]);
    }
    //=======================================================================================================

    int nT=RealSourcesData.ncol();
    matrice* data;

    if(!strcmp(argv[1],"-EEG"))
    {
        SimulatedEegData=matrice(EegGainMatrix.nlin(),nT);
        data=&SimulatedEegData;
    }
    if(!strcmp(argv[1],"-MEG"))
    {
        SimulatedMegData=matrice(MegGainMatrix.nlin(),nT);
        data=&SimulatedMegData;
    }

    double noiselevel;
    
    #ifdef USE_OMP
    #pragma omp parallel for
    #endif
    for(int frame=0;frame<nT;frame++)
    {
        vecteur v; v.DangerousBuild(&RealSourcesData(0,frame),RealSourcesData.nlin());
        vecteur result;

        if(!strcmp(argv[1],"-EEG"))
        {
            result=EegGainMatrix*v;
            noiselevel=EegNoiseLevel;
        }
        if(!strcmp(argv[1],"-MEG"))
        {
            result=MegGainMatrix*v;
            noiselevel=MegNoiseLevel;
        }
        for(size_t i=0;i<result.size();i++) (*data)(i,frame)=result(i);
        v.DangerousKill();
    }

    double meannorm = 0.0;

    for(int frame=0;frame<nT;frame++)
    {
        vecteur v; v.DangerousBuild(&((*data)(0,frame)),data->nlin());
        meannorm += (1.0/nT)*v.norm();
        v.DangerousKill();
    }

    for(int frame=0;frame<nT;frame++)
    {
        vecteur grand(data->nlin());
        for(size_t i=0;i<grand.size();i++) grand(i)=gaussienne();
        double coef=noiselevel*meannorm/grand.norm();

        vecteur v; v.DangerousBuild(&((*data)(0,frame)),data->nlin());
        for(size_t i=0;i<v.size();i++) v(i)+=coef*grand(i);
        v.DangerousKill();
    }

    // write output variables ===================================================================================
    // for use with EEG DATA
    if(!strcmp(argv[1],"-EEG"))
    {
        SimulatedEegData.saveTxt(argv[4]);
    }
    // for use with MEG DATA
    else if(!strcmp(argv[1],"-MEG"))
    {
        SimulatedMegData.saveTxt(argv[4]);
    }
    // ===========================================================================================================

    // Stop Chrono
    C.stop();
    C.dispEllapsed();

    return 0;
}

void getHelp(char** argv)
{
    cout << argv[0] <<" [-option] [filepaths...]" << endl << endl;
    
    cout << "-option :" << endl;
    cout << "   -EEG :   Compute the forward problem for EEG " << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            EegGainMatrix (bin), RealSourcesData (txt), SimulatedEegData (txt), EegNoiseLevel (float)" << endl << endl;

    cout << "   -MEG :   Compute the forward problem for MEG " << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            MegGainMatrix (bin), RealSourcesData (txt), SimulatedMegData (txt), MegNoiseLevel (float)" << endl << endl;

    exit(0);
}


