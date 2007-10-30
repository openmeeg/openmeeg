#include "symmatrice.h"
#include "vecteur.h"
#include "mesh3.h"
#include "sparse_matrice.h"
#include "fast_sparse_matrice.h"

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

    // declaration of argument variables
    Mesh SourceMesh;

    bool checkClosedSurface = false;
    SourceMesh.load(argv[1],checkClosedSurface);

    sparse_matrice SmoothMatrix = SourceMesh.gradient();
    vecteur AiVector = SourceMesh.areas();

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
    cout << "   SourceMesh, SmoothMatrix (bin), AiVector (bin)" << endl << endl;

    exit(0);
}

