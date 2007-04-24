#include "MatLibConfig.h"
#include "vecteur.h"
#include "matrice.h"
#include "symmatrice.h"
#include "sparse_matrice.h"
#include "fast_sparse_matrice.h"
#include <iostream>
#include <cmath>

using namespace std;

template<class T> bool compare(char** argv, T mat1, T mat2, float eps=0.00001);
template<class T> double normInf(T mat);
void getHelp(char** argv);

int main (int argc, char** argv)
{
    if(argc < 2)
    {
        cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
        return 0;
    }
    if ((!strcmp(argv[1],"-h")) | (!strcmp(argv[1],"--help"))) getHelp(argv);

    if(argc < 5)
    {
        cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
        return 0;
    }
    float eps=0.0001; // arbitrary
    bool flag;
    if(strcmp(argv[2],"-sym")==0){
        symmatrice M,Q;
        flag=compare<symmatrice>(argv,M,Q,eps);
    }
    else
        if(strcmp(argv[2],"-std")==0){
            matrice M,Q;
            flag=compare<matrice>(argv,M,Q,eps);
        }
        else{
            cerr << "Wrong argument ! \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
            exit(1);
        }

    if(!flag){
        cerr << endl << "ERROR : Matrices are different at the precision : eps = " << eps << endl;
        exit(1);
    }
    cout << "OK" << endl;
    return 0;
}


template<class T>
bool compare(char** argv, T mat1, T mat2, float eps){
// T is a matrice or a symmatrice

    cout << "-------------------------------------------------" << endl;
    cout << "Comparing : " << endl;
    cout << "- " << argv[3] << endl;
    cout << "- " << argv[4] << endl;

    if(strcmp(argv[1],"-bin")==0) {
        mat1.loadBin(argv[3]);
        mat2.loadBin(argv[4]);
    }
    else if(strcmp(argv[1],"-txt")==0) {
        mat1.loadTxt(argv[3]);
        mat2.loadTxt(argv[4]);
    }
    else {
        cerr << "Wrong argument ! \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
        exit(1);
    }

    if ((mat1.ncol() != mat2.ncol()) || (mat1.nlin() != mat2.nlin())) {
        cerr << "ERROR : Dimension mismatch !" << endl;
        exit(1);
    }

    bool flag=true;

    double norm1 = normInf<T>(mat1);
    double norm2 = normInf<T>(mat2);
    double diff;
    if((norm1>1e-4)&(norm2>1e-4)) {
        for(unsigned int i=0; i<mat1.nlin(); i++) {
            for(unsigned int j=0; j<mat1.ncol(); j++) {
                diff = abs(mat1(i,j)/norm1 - mat2(i,j)/norm2);
                flag = flag && (diff < eps);
                if(!(abs(mat1(i,j)/norm1 - mat2(i,j)/norm2) < eps)) {
                    cout << "NORM  " << mat1(i,j) << "  " << mat2(i,j) << "  " << diff << endl;
                }
            }
        }
    } else {
        for(unsigned int i=0; i<mat1.nlin(); i++) {
            for(unsigned int j=0; j<mat1.ncol(); j++) {
                if(abs(mat1(i,j))>1e-4){
                    flag = flag && (abs(mat1(i,j) - mat2(i,j))/abs(mat1(i,j)) < eps);
                    if(!(abs(mat1(i,j) - mat2(i,j))/mat1(i,j) < eps))
                        cout << " ERR REL  " << mat1(i,j) << "  " << mat2(i,j) << "  " << abs(mat1(i,j) - mat2(i,j))/mat1(i,j) << endl;
                }
                else
                    flag = flag && (abs(mat1(i,j) - mat2(i,j)) < eps);
            }
        }
    }
    return flag;
}

template<class T>
double normInf(T mat){ // compute the max of the norm 1 of each line
    double max = 0.0;
    double sum;
    for(unsigned int i=0;i<mat.nlin();i++){
        sum = 0.0;
        for(unsigned int j=0;j<mat.ncol();j++) {
            sum += abs(mat(i,j));
        }
        if(max < sum)
            max=sum;
    }
    return max;
}

void getHelp(char** argv)
{
    cout << argv[0] <<" [-options : -(bin|txt) -(sym|std)] [filepaths...]" << endl << endl;

    cout << "-option :" << endl;
    cout << "   -bin :   Compare matrices written in binary format" << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            matrix1 , matrix2  " << endl << endl;

    cout << "   -txt :   Compare matrices written in text format" << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            matrix1 , matrix2  " << endl << endl;

    cout << "   -sym :   Compare symmetric matrices" << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            matrix1 , matrix2  " << endl << endl;

    cout << "   -std :   Compare standard matrices" << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            matrix1 , matrix2  " << endl << endl;

    exit(0);
}
