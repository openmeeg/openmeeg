#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <matrice.h>
#include <symmatrice.h>
#include <vecteur.h>

using namespace std;
//using namespace CLMatLib;

int main( int argc, char **argv)
{
    if(argc!=3 && argc!=4 || argc==4 && strcmp(argv[3],"--abs") )
    {
        cout << " usage: "<< argv[0] <<" input_txt_file output_texture_file [--abs]" << endl;
        exit(1);
    }

	matrice M(argv[1],'t');
    ofstream os(argv[2]);
    os<<"ascii"<<endl<<"FLOAT"<<endl<<(unsigned int)M.ncol()<<endl;
    for(size_t j=0;j<M.ncol();j++)
    {
        os<<(unsigned int)j<<endl;
        os<<(unsigned int)M.nlin();
        if(argc==3) for(size_t i=0;i<M.nlin();i++) os<<" "<<M(i,j);
        else for(size_t i=0;i<M.nlin();i++) os<<" "<<fabs(M(i,j));
        os<<endl;
    }

    os.close();

    return 0;
}
