#include <matrice.h>
#include <symmatrice.h>
#include <vecteur.h>

using namespace std;
//using namespace CLMatLib;

int main( int argc, char **argv)
{
    if(argc!=3)
    {
        cout << "usage: " << argv[0] << " input_matrix_txt_file output_matrix_bin_file" << endl;
        exit(1);
    }

    symmatrice M(argv[1],'t');
    M.saveBin(argv[2]);
    return 0;
}
