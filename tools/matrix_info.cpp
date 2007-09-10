#include <vecteur.h>
#include <symmatrice.h>
#include <matrice.h>
#include <cmath>

#include "options.h"

using namespace std;

template<class T> void print_infos(const T& M);

int main( int argc, char **argv)
{
    command_usage("Compare two matrices of float with a certain numerical precision\ncompare_matrix mat1 mat2 [options]");
    const char *filename = command_option("-i",(const char *) NULL,"Matrix file");
    const char *txt = command_option("-txt",(const char *) 0,"Force reading data stored in ascii format");
    const char *sym = command_option("-sym",(const char *) 0,"Data are symmetric matrices");
    if (command_option("-h",(const char *)0,0)) return 0;
    
    if(!filename)
    {
        cerr << "Please set Matrix File" << endl;
        exit(1);
    }

    cout << "Loading : " << filename << endl;

    if(sym) {
        if(txt)
        {
            symmatrice M(filename,'t');
            cout << "Format : ASCII" << endl;
            print_infos(M);
        } else {
            symmatrice M(filename,'b');
            cout << "Format : BINARY" << endl;
            print_infos(M);
        }
    } else {
        if(txt)
        {
            matrice M(filename,'t');
            cout << "Format : ASCII" << endl;
            print_infos(M);
        } else {
            matrice M(filename,'b');
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
    for(size_t i = 0; i < std::min(M.nlin(),(size_t) 5); ++i)
    {
        for(size_t j = 0; j < std::min(M.ncol(),(size_t) 5); ++j)
        {
            cout << M(i,j) << " " ;
        }
        cout << endl ;
    }
}
