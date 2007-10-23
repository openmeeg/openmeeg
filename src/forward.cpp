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

    cout << endl << "| ------ " << argv[0] << " -------" << endl;
    for( int i = 1; i < argc; i += 1 )
    {
        cout << "| " << argv[i] << endl;
    }
    cout << "| -----------------------" << endl;

    // declaration of argument variables======================================================================
    string Option;
    matrice GainMatrix;
    matrice RealSourcesData;
    double NoiseLevel;
    matrice SimulatedData;

    GainMatrix.loadBin(argv[1]);
    RealSourcesData.loadTxt(argv[2]);
    NoiseLevel = atof(argv[4]);

    int nT = RealSourcesData.ncol();
    matrice* data;

    SimulatedData=matrice(GainMatrix.nlin(),nT);
    data = &SimulatedData;

    double noiselevel;

    #ifdef USE_OMP
    #pragma omp parallel for
    #endif
    for(int frame=0;frame<nT;frame++)
    {
        vecteur v; v.DangerousBuild(&RealSourcesData(0,frame),RealSourcesData.nlin());
        vecteur result;

        result = GainMatrix*v;
        noiselevel = NoiseLevel;

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
        for(size_t i=0;i<grand.size();i++) grand(i)=gaussian();
        double coef=noiselevel*meannorm/grand.norm();

        vecteur v; v.DangerousBuild(&((*data)(0,frame)),data->nlin());
        for(size_t i=0;i<v.size();i++) v(i)+=coef*grand(i);
        v.DangerousKill();
    }

    // write output variables ===================================================================================
    SimulatedData.saveTxt(argv[3]);
    // ===========================================================================================================

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


