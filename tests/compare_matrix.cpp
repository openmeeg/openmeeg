#include "MatLibConfig.h"
#include "vector.h"
#include "matrix.h"
#include "symmatrix.h"
#include "sparse_matrix.h"
#include "options.h"
#include <iostream>
#include <cmath>

using namespace std;

template<class T> bool compare(const T& mat1, const T& mat2, double eps, size_t col = 0);
template<class T> bool compare_rdm(const T& mat1, const T& mat2, double eps, size_t col = 0);
template<class T> double normInf(const T& mat);
template<class T> bool compare_matrixs(maths::ifstream& ifs1,T& mat1,maths::ifstream& ifs2,T& mat2,
                                          double eps,const char* rdm,size_t col);

int main (int argc, char** argv)
{
    command_usage("Compare two matrices with a certain numerical precision\ncompare_matrix mat1 mat2 [options]");
    const char *input_format1 = command_option("-if1",(const char *) NULL,
                                                "Input file format for Matrix 1 : ascii, binary, tex, matlab");
    const char *input_format2 = command_option("-if2",(const char *) NULL,
                                                "Input file format for Matrix 2 : ascii, binary, tex, matlab");

    const char *isfull = command_option("-full",(const char *) 0,"Data are symmetric matrices");
    const char *issym = command_option("-sym",(const char *) 0,"Data are symmetric matrices");
    const char *issparse = command_option("-sparse",(const char *) 0,"Data are sparse matrices");
    const char *epsilon = command_option("-eps","0.00001","Tolerance on differences");
    const char *rdm = command_option("-rdm",(const char *) 0,"Use RDM (Relative difference measure) to compare each column of matrices");
    const int col = command_option("-col",(int) 0,"Restrict RDM comparison to one column (index starts at 1)");
    if (command_option("-h",(const char *)0,0)) return 0;

    if(argc < 3) {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    if(!isfull && !issym && !issparse) {
        std::cout << "Please set Matrix type using : -full, -sym or -sparse" << std::endl;
        return 1;
    }

    if(!isfull && rdm) {
        std::cerr << "Can use -rdm or -col only with full matrices" << std::endl;
        return 1;
    }

    double eps = atof(epsilon);

    cout << "-------------------------------------------------" << endl;
    cout << "Comparing : " << endl;
    cout << "- " << argv[1] << endl;
    cout << "- " << argv[2] << endl;

    maths::ifstream ifs1(argv[1]);
    maths::ifstream ifs2(argv[2]);

    if(input_format1) {
        ifs1 = ifs1 >> maths::format(input_format1);
    }

    if(input_format2) {
        ifs2 = ifs2 >> maths::format(input_format2);
    }

    bool flag;
    if(issym) {
        SymMatrix mat1;
        SymMatrix mat2;
        ifs1 >> mat1;
        ifs2 >> mat2;
        flag = compare(mat1,mat2,eps,col);
    } else if (issparse) {
        SparseMatrix mat1;
        SparseMatrix mat2;
        ifs1 >> mat1;
        ifs2 >> mat2;
        flag = compare(mat1,mat2,eps,col);
    } else { // assumes isfull
        Matrix mat1;
        Matrix mat2;
        flag = compare_matrixs(ifs1,mat1,ifs2,mat2,eps,rdm,col);
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
bool compare_matrixs(maths::ifstream& ifs1,T& mat1,maths::ifstream& ifs2,T& mat2,
                                          double eps,const char*rdm,size_t col) {
    bool flag;
    try
    {
        ifs1 >> mat1;
        ifs2 >> mat2;
    } catch (std::string s) {
        std::cerr << s << std::endl;
        exit(1);
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
    return flag;
}

template<class T>
bool compare(const T& mat1, const T& mat2, double eps, size_t col){
// T is a Matrix or a SymMatrix

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
bool compare_rdm(const T& mat1, const T& mat2, double eps, size_t col){
// T is a Matrix

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
        Vector col1 = mat1.getcol(j);
        Vector col2 = mat2.getcol(j);
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
