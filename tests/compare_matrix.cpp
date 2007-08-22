#include "MatLibConfig.h"
#include "vecteur.h"
#include "matrice.h"
#include "symmatrice.h"
#include "sparse_matrice.h"
#include "fast_sparse_matrice.h"
#include "options.h"
#include <iostream>
#include <cmath>

using namespace std;

template<class T> bool compare(const T& mat1, const T& mat2, float eps);
template<class T> bool compare_corr(const T& mat1, const T& mat2, float eps);
template<class T> double normInf(const T& mat);

int main (int argc, char** argv)
{
    if (argc<3) {
        cerr << "Error: Not enough arguments, try the -h option" << endl;
        return 1;
    }

    command_usage("Compare two matrices of float with a certain numerical precision\ncompare_matrix mat1 mat2 [options]");
    const char *bin = command_option("-bin",(const char *) 0,"Force reading data stored in binary format");
    const char *txt = command_option("-txt",(const char *) 0,"Force reading data stored in ascii format");
    const char *sym = command_option("-sym",(const char *) 0,"Data are symmetric matrices");
    const char *epsilon = command_option("-eps","0.0001","Data are symmetric matrices");
    const char *corr = command_option("-corr",(const char *) 0,"Use 1-correlation between normalized columns of matrices");
    if (command_option("-h",(const char *)0,0)) return 0;

    float eps = atof(epsilon);

    cout << "-------------------------------------------------" << endl;
    cout << "Comparing : " << endl;
    cout << "- " << argv[1] << endl;
    cout << "- " << argv[2] << endl;
    
    bool flag;
    if(sym){
        symmatrice mat1;
        symmatrice mat2;
        if(bin) {
            mat1.loadBin(argv[1]);
            mat2.loadBin(argv[2]);
        } else if(txt) {
            mat1.loadTxt(argv[1]);
            mat2.loadTxt(argv[2]);
        } else {
            mat1.load(argv[1]);
            mat2.load(argv[2]);
        }
        if(corr) {
            std::cerr << "ERROR : Cannot use correlation on symmetric matrices" << std::endl;
            exit(1);
        } else {
            flag = compare(mat1,mat2,eps);
        }
    }
    else {
        matrice mat1;
        matrice mat2;
        if(bin) {
            mat1.loadBin(argv[1]);
            mat2.loadBin(argv[2]);
        } else if(txt) {
            mat1.loadTxt(argv[1]);
            mat2.loadTxt(argv[2]);
        } else {
            mat1.load(argv[1]);
            mat2.load(argv[2]);
        }
        if(corr) {
            flag = compare_corr(mat1,mat2,eps);
        } else {
            flag = compare(mat1,mat2,eps);
        }
    }

    if(!flag){
        cerr << endl << "ERROR : Matrices are different at the precision : eps = " << eps << endl;
        exit(1);
    }
    cout << "OK" << endl;
    return 0;
}

template<class T>
bool compare(const T& mat1, const T& mat2, float eps){
// T is a matrice or a symmatrice

    if ((mat1.ncol() != mat2.ncol()) || (mat1.nlin() != mat2.nlin())) {
        cerr << "ERROR : Dimension mismatch !" << endl;
        exit(1);
    }

    bool flag = true;

    double norm1 = normInf<T>(mat1);
    double norm2 = normInf<T>(mat2);
    double diff;

    if((norm1>1e-4)&(norm2>1e-4)) {
        for(unsigned int i=0; i<mat1.nlin(); i++) {
            for(unsigned int j=0; j<mat1.ncol(); j++) {
                diff = abs(mat1(i,j)/norm1 - mat2(i,j)/norm2);
                flag = flag && (diff < eps);
                if(!(diff < eps)) {
                    cout << "NORM  " << mat1(i,j) << "  " << mat2(i,j) << "  " << diff << endl;
                }
            }
        }
    } else {
        for(unsigned int i=0; i<mat1.nlin(); i++) {
            for(unsigned int j=0; j<mat1.ncol(); j++) {
                if(abs(mat1(i,j))>1e-4){
                    diff = abs(mat1(i,j) - mat2(i,j))/abs(mat1(i,j));
                    flag = flag && ( diff < eps);
                    if(!(diff < eps))
                        cout << " ERR REL  " << mat1(i,j) << "  " << mat2(i,j) << "  " << diff << endl;
                }
                else
                    flag = flag && (abs(mat1(i,j) - mat2(i,j)) < eps);
            }
        }
    }
    return flag;
}

template<class T>
bool compare_corr(const T& mat1, const T& mat2, float eps){
// T is a matrice

    if ((mat1.ncol() != mat2.ncol()) || (mat1.nlin() != mat2.nlin())) {
        cerr << "ERROR : Dimension mismatch !" << endl;
        exit(1);
    }

    bool flag = true;
    double diff;
    for(unsigned int j=0; j<mat1.ncol(); j++) {
        vecteur col1 = mat1.getcol(j);
        vecteur col2 = mat2.getcol(j);
        col1 = col1 - col1.mean();
        col2 = col2 - col2.mean();
        col1 = col1 / col1.norm();
        col2 = col2 / col2.norm();
        diff = (col1 - col2).norm() / mat1.nlin();
        flag = flag && (diff < eps);
        if(!(diff < eps)) {
            cout << "NORM ( column " << j << " ) " << diff << endl;
        }
    }
    return flag;
}

template<class T>
double normInf(const T& mat){ // compute the max of the norm 1 of each line
    double max = 0.0;
    double sum;
    for(unsigned int i=0;i<mat.nlin();i++){
        sum = 0.0;
        for(unsigned int j=0;j<mat.ncol();j++) {
            sum += abs(mat(i,j));
        }
        if(max < sum)
            max = sum;
    }
    return max;
}
