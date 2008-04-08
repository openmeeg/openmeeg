#include "MatLibConfig.h"
#include "vecteur.h"
#include "matrice.h"
#include "symmatrice.h"
#include "sparse_matrice.h"
#include "options.h"
#include <iostream>
#include <cmath>

using namespace std;

template<class T> bool compare(const T& mat1, const T& mat2, float eps, size_t col = 0);
template<class T> bool compare_rdm(const T& mat1, const T& mat2, float eps, size_t col = 0);
template<class T> double normInf(const T& mat);

int main (int argc, char** argv)
{
    command_usage("Compare two matrices of float with a certain numerical precision\ncompare_matrix mat1 mat2 [options]");
    const char *bin = command_option("-bin",(const char *) 0,"Force reading data stored in binary format");
    const char *txt = command_option("-txt",(const char *) 0,"Force reading data stored in ascii format");
    const char *sym = command_option("-sym",(const char *) 0,"Data are symmetric matrices");
    const char *sparse = command_option("-sparse",(const char *) 0,"Data are sparse matrices");
    const char *epsilon = command_option("-eps","0.00001","Tolerance on differences");
    const char *rdm = command_option("-rdm",(const char *) 0,"Use RDM (Relative difference measure) to compare each column of matrices");
    const int col = command_option("-col",(int) 0,"Restrict RDM comparison to one column (index starts at 1)");
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
        if(rdm) {
            std::cerr << "ERROR : Cannot use RDM on symmetric matrices" << std::endl;
            exit(1);
        } else {
            flag = compare(mat1,mat2,eps);
        }
    } else if (sparse) {
        sparse_matrice mat1;
        sparse_matrice mat2;
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
        if(rdm) {
            std::cerr << "ERROR : Cannot use RDM on sparse matrices" << std::endl;
            exit(1);
        } else {
            flag = compare(mat1,mat2,eps);
        }
    } else {
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
        if(rdm) {
            if(col) {
                flag = compare_rdm(mat1,mat2,eps,col);
            } else {
                flag = compare_rdm(mat1,mat2,eps);
            }
        } else {
            if(col) {
                flag = compare(mat1,mat2,eps,col);
            } else {
                flag = compare(mat1,mat2,eps);
            }
        }
    }

    if(!flag){
        cerr << endl << "ERROR : Matrices are different at the precision : eps = " << eps << endl;
        exit(1);
    }
    cout << "OK" << endl;
    cout.flush();
    return 0;
}

template<class T>
bool compare(const T& mat1, const T& mat2, float eps, size_t col){
// T is a matrice or a symmatrice

    if(col) {
        if ((mat1.ncol() < col) || (mat2.ncol() < col)) {
            cerr << "ERROR : Bad Column Id for matrices dimensions !" << endl;
            exit(1);
        }
    } else {
        if ((mat1.ncol() != mat2.ncol()) || (mat1.nlin() != mat2.nlin())) {
            cerr << "ERROR : Dimension mismatch !" << endl;
            exit(1);
        }
    }

    unsigned int jmin,jmax;
    if(col > 0) {
        jmin = col-1;
        jmax = col;
    } else {
        jmin = 0;
        jmax = mat1.ncol();;
    }

    bool flag = true;

    double norm1 = normInf<T>(mat1);
    double norm2 = normInf<T>(mat2);
    double diff;

    if((norm1>1e-4)&(norm2>1e-4)) {
        for(unsigned int i=0; i<mat1.nlin(); i++) {
            for(unsigned int j=jmin; j<jmax; j++) {
                diff = abs(mat1(i,j)/norm1 - mat2(i,j)/norm2);
                flag = flag && (diff < eps);
                if (!(diff < eps)) {
                    cout << "ERROR NORM  " << mat1(i,j) << "  " << mat2(i,j) << "  " << diff << endl;
                    cout.flush();
                }
            }
        }
    } else {
        for(unsigned int i=0; i<mat1.nlin(); i++) {
            for(unsigned int j=jmin; j<jmax; j++) {
                if (abs(mat2(i,j))>1e-4) {
                    diff = abs(mat1(i,j) - mat2(i,j))/abs(mat2(i,j));
                    flag = flag && ( diff < eps);
                    if (!(diff < eps))
                        cout << "ERROR RELATIVE  " << mat1(i,j) << "  " << mat2(i,j) << "  " << diff << endl;
                        cout.flush();
                }
                else {
                    diff = abs(mat1(i,j) - mat2(i,j));
                    flag = flag && ( diff < eps);
                    if (!(diff < eps)) {
                        cout << "ERROR DIFF  " << mat1(i,j) << "  " << mat2(i,j) << "  " << diff << endl;
                        cout.flush();
                    }
                }
            }
        }
    }
    return flag;
}

template<class T>
bool compare_rdm(const T& mat1, const T& mat2, float eps, size_t col){
// T is a matrice

    if(col) {
        if ((mat1.ncol() < col) || (mat2.ncol() < col)) {
            cerr << "ERROR : Bad Column Id for matrices dimensions !" << endl;
            exit(1);
        }
    } else {
        if ((mat1.ncol() != mat2.ncol()) || (mat1.nlin() != mat2.nlin())) {
            cerr << "ERROR : Dimension mismatch !" << endl;
            exit(1);
        }
    }

    bool flag = true;
    double diff;
    unsigned int jmin,jmax;
    if(col > 0) {
        jmin = col-1;
        jmax = col;
    } else {
        jmin = 0;
        jmax = mat1.ncol();;
    }

    for(unsigned int j=jmin; j<jmax; j++) {
        vecteur col1 = mat1.getcol(j);
        vecteur col2 = mat2.getcol(j);
        col1 = col1 - col1.mean();
        col2 = col2 - col2.mean();
        col1 = col1 / col1.norm();
        col2 = col2 / col2.norm();
        diff = (col1 - col2).norm();

        flag = flag && (diff < eps);
        if(diff > eps) {
            cout << "ERROR RDM ( column " << j << " ) " << diff << endl;
            cout.flush();
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
