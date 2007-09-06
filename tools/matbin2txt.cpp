#include <vecteur.h>
#include <symmatrice.h>
#include <matrice.h>

using namespace std;

int main( int argc, char **argv)
{
    if(argc!=3)
    {
        cout << "usage: " << argv[0] << " input_matrix_file output_matrix_file" << endl;
        exit(1);
    }

    matrice M(argv[1],'b');
    M.saveTxt(argv[2]);
    
    return 0;
}
