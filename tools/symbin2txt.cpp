/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
*/

#include <vecteur.h>
#include <symmatrice.h>
#include <matrice.h>

using namespace std;

int main( int argc, char **argv)
{
    if(argc!=3)
    {
        cout << "usage: " << argv[0] << " input_matrix_bin_file output_matrix_txt_file" << endl;
        exit(1);
    }

    symmatrice M(argv[1],'b');
    M.saveTxt(argv[2]);
    return 0;
}
