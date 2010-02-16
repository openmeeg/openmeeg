/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Maureen.Clerc.AT.sophia.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.sophia.inria.fr)

The OpenMEEG software is a C++ package for solving the forward/inverse
problems of electroencephalography and magnetoencephalography.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's authors,  the holders of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

#include "vector.h"
#include "sparse_matrix.h"
#include "symmatrix.h"
#include "matrix.h"
#include "options.h"

#include <cmath>

using namespace std;
using namespace OpenMEEG;

template<class T> void print_infos(const T& M);

int main( int argc, char **argv)
{
    print_version(argv[0]);

    command_usage("Provides informations on a Matrix generated with OpenMEEG");
    const char *filename = command_option("-i",(const char *) NULL,"Matrix file");
    const char *txt = command_option("-txt",(const char *) 0,"Force reading data stored in ascii format");
    const char *sym = command_option("-sym",(const char *) 0,"Data are symmetric matrices");
    const char *sparse = command_option("-sparse",(const char *) 0,"Data are sparse matrices");
    const char *mat = command_option("-mat",(const char *) 0,"Data are matlab format");
    if (command_option("-h",(const char *)0,0)) return 0;

    if(!filename)
    {
        cerr << "Please set Matrix File" << endl;
        exit(1);
    }

    cout << "Loading : " << filename << endl;

    if(sym) {
        if(txt) {
            SymMatrix M;
            M.load(filename);
            cout << "Format : ASCII" << endl;
            print_infos(M);
        } else if(mat) {
            cerr << "Unsupported Format : MAT for symmetric matrices" << endl;
            exit(1);
        } else {
            SymMatrix M;
            M.load(filename);
            cout << "Format : BINARY" << endl;
            print_infos(M);
        }
    } else if(sparse) {
        if(txt) {
            SparseMatrix M;
            M.load(filename);
            cout << "Format : ASCII" << endl;
            print_infos(M);
        } else if(mat) {
            SparseMatrix M;
            M.load(filename);
            cout << "Format : MAT" << endl;
            print_infos(M);
        } else {
            SparseMatrix M;
            M.load(filename);
            cout << "Format : BINARY" << endl;
            print_infos(M);
        }
    } else {
        if(txt) {
            Matrix M;
            M.load(filename);
            cout << "Format : ASCII" << endl;
            print_infos(M);
        } else if(mat) {
            Matrix M;
            M.load(filename);
            cout << "Format : MAT" << endl;
            print_infos(M);
        } else {
            Matrix M;
            M.load(filename);
            cout << "Format : BINARY" << endl;
            print_infos(M);
        }
    }

    return 0;
}

template<class T>
void print_infos(const T& M) {
    if ((M.nlin() == 0) && (M.ncol() == 0)) {
        cout << "Matrix Empty" << endl;
        return;
    }

    cout << "Dimensions : " << M.nlin() << " x " << M.ncol() << endl;

    double minv = M(0,0);
    double maxv = M(0,0);
    size_t mini = 0;
    size_t maxi = 0;
    size_t minj = 0;
    size_t maxj = 0;

    for(size_t i = 0; i < M.nlin(); ++i)
    {
        for(size_t j = 0; j < M.ncol(); ++j)
        {
            if (minv > M(i,j)) {
                minv = M(i,j);
                mini = i;
                minj = j;
            } else if (maxv < M(i,j)) {
                maxv = M(i,j);
                maxi = i;
                maxj = j;
            }
        }
    }
    cout << "Min Value : " << minv << " (" << mini << "," << minj << ")" << endl;
    cout << "Max Value : " << maxv << " (" << maxi << "," << maxj << ")" << endl;

    cout << "First Values" << endl;
    for(size_t i = 0; i < std::min(M.nlin(),(size_t) 10); ++i)
    {
        for(size_t j = 0; j < std::min(M.ncol(),(size_t) 10); ++j)
        {
            cout << M(i,j) << " " ;
        }
        cout << endl ;
    }
}
