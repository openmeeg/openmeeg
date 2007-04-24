#include "symmatrice.h"
#include "vecteur.h"
#include "sparse_matrice.h"
#include "fast_sparse_matrice.h"

#include "MeshDataSmoother.h"

using namespace std;
//using namespace CLMatLib;

void getHelp(char** argv);

int main(int argc, char **argv)
{
    if(argc==1)
    {
        cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
        return 0;
    }

    if ((!strcmp(argv[1],"-h")) | (!strcmp(argv[1],"--help"))) getHelp(argv);

    // declaration of argument variables
    mesh<3> SourceMesh;
    sparse_matrice SmoothMatrix;
    fast_sparse_matrice fastSmoothMatrix;
    fast_sparse_matrice fastSmoothMatrix_t;
    vecteur AiVector;
    
    SourceMesh.load(argv[1]);

    MeshDataL1Phi MD2(SourceMesh);
    std::vector<double> *Ai;
    MD2.computeMatrix(Ai);
    sparse_matrice *mat=MD2.getMatrix();
    mat->refreshNZ();
    SmoothMatrix=(*mat);
    vecteur v;
    v.DangerousBuild(&(*Ai)[0],Ai->size());
    AiVector=v.duplicate();
    v.DangerousKill();
    delete Ai;

    // write output variables
    SmoothMatrix.saveBin(argv[2]);
    AiVector.saveBin(argv[3]);

    return 0;
}

void getHelp(char** argv)
{
    cout << argv[0] <<" [filepaths...]" << endl << endl;

    cout << "   Presmooth a mesh " << endl;
    cout << "   Filepaths are in order :" << endl;
    cout << "   SourceMesh, SmoothMatrix (txt), AiVector (txt)" << endl << endl;

    exit(0);

}


