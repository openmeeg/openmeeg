#include "OpenMEEGMathsConfig.h"
#include "vector.h"
#include "matrix.h"
#include "symmatrix.h"
#include "sparse_matrix.h"
#include "options.h"
#include <iostream>
#include <cmath>

using namespace OpenMEEG;

template<class T> bool compare(const T& mat1, const T& mat2, double eps, size_t col = 0);
template<class T> bool compare_rdm(const T& mat1, const T& mat2, double eps, size_t col = 0);
template<class T> bool compare_mag(const T& mat1, const T& mat2, double eps, size_t col = 0);
template<class T> bool compare_matrix(maths::ifstream& ifs1,T& mat1,maths::ifstream& ifs2,T& mat2,
        double eps,const char* rdm,const char* mag,size_t col);
bool compare(maths::ifstream& ifs1,Vector& V1,maths::ifstream& ifs2,Vector& V2,const double eps);

int main (int argc, char** argv)
{
    command_usage("Compare two matrices with a certain numerical precision\ncompare_matrix mat1 mat2 [options]");
    // const char *input_format1 = command_option("-if1",(const char *) NULL,
    //         "Input file format for Matrix 1 : ascii, binary, tex, matlab");
    // const char *input_format2 = command_option("-if2",(const char *) NULL,
    //         "Input file format for Matrix 2 : ascii, binary, tex, matlab");

    const char* isfull   = command_option("-full",(const char *) 0,"Data are full matrices");
    const char* issym    = command_option("-sym",(const char *) 0,"Data are symmetric matrices");
    const char* issparse = command_option("-sparse",(const char *) 0,"Data are sparse matrices");
    const char* isvector = command_option("-vector",(const char *) 0,"Data are vectors");

    const char* epsilon = command_option("-eps","0.00002","Tolerance on differences"); // Hacking tol for tests

    const char* rdm = command_option("-rdm",(const char *) 0,"Use RDM (Relative difference measure) to compare each column of matrices");
    const char* mag = command_option("-mag",(const char *) 0,"Use MAG (MAGnification error) to compare each column of matrices");

    const int col = command_option("-col",(int) 0,"Restrict RDM comparison to one column (index starts at 1)");

    if (command_option("-h",(const char *)0,0)) return 0;

    if (argc<3) {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    if (!isfull && !issym && !issparse && !isvector) {
        std::cout << "Please set Matrix type using : -vector, -full, -sym or -sparse" << std::endl;
        return 1;
    }

    if (!isfull && (rdm||mag)) {
        std::cerr << "Can use -rdm, -mag or -col only with full matrices" << std::endl;
        return 1;
    }

    if (rdm && mag) {
        std::cerr << "Choose either -rdm OR -mag but not both" << std::endl;
        return 1;
    }


    double eps = atof(epsilon);

    std::cout << "-------------------------------------------------" << std::endl;
    std::cout << "Comparing : " << std::endl;
    std::cout << "- " << argv[1] << std::endl;
    std::cout << "- " << argv[2] << std::endl;

    maths::ifstream ifs1(argv[1]);
    if (!ifs1) {
        std::cerr << "Cannot open file " << argv[1] << " !" << std::endl;
        return 1;
    }

    maths::ifstream ifs2(argv[2]);
    if (!ifs2) {
        std::cerr << "Cannot open file " << argv[2] << " !" << std::endl;
        return 1;
    }

#if 0
    if (input_format1) {
        ifs1 >> maths::format(input_format1);
    }

    if (input_format2) {
        ifs2 >> maths::format(input_format2);
    }
#endif

    bool flag;
    if (issym) {
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
    } else if (isfull) {
        Matrix mat1;
        Matrix mat2;
        flag = compare_matrix(ifs1,mat1,ifs2,mat2,eps,rdm,mag,col);
    } else {
        Vector V1;
        Vector V2;
        flag = compare(ifs1,V1,ifs2,V2,eps);
    }

    if (!flag){
        std::cerr << std::endl << "ERROR : Matrices are different at the precision : eps = " << eps << std::endl;
        exit(1);
    }
    std::cout << "OK" << std::endl;
    std::cout.flush();

    return 0;
}

template<class T>
double normInf(const T& mat){ // compute the max of the norm 1 of each line
    double max = 0.0;
    double sum;
    for(unsigned int i=0;i<mat.nlin();++i){
        sum = 0.0;
        for(unsigned int j=0;j<mat.ncol();++j) {
            sum += std::abs(mat(i,j));
        }
        if (max<sum)
            max = sum;
    }
    return max;
}

double normInf(const Vector& V) { // compute the norm 1 of a vector
    double sum = 0.0;
    for (unsigned i=0;i<V.size();++i)
        sum += std::abs(V(i));
    return sum;
}

template<class T>
bool compare_matrix(maths::ifstream& ifs1,T& mat1,maths::ifstream& ifs2,T& mat2,
        double eps,const char*rdm,const char*mag,size_t col) {
    bool flag;
    try {
        ifs1 >> mat1;
        ifs2 >> mat2;
    } catch (std::string s) {
        std::cerr << s << std::endl;
        exit(1);
    }
    
    if (rdm) {
        if (col) {
            flag = compare_rdm(mat1,mat2,eps,col);
        } else {
            flag = compare_rdm(mat1,mat2,eps);
        }
    } else if (mag) {
        if (col) {
            flag = compare_mag(mat1,mat2,eps,col);
        } else {
            flag = compare_mag(mat1,mat2,eps);
        }
    } else {
        if (col) {
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

    if (col) {
        if ((mat1.ncol()<col) || (mat2.ncol()<col)) {
            std::cerr << "ERROR : Bad Column Id for matrices dimensions !" << std::endl;
            exit(1);
        }
    } else {
        if ((mat1.ncol()!=mat2.ncol()) || (mat1.nlin()!=mat2.nlin())) {
            std::cerr << "ERROR : Dimension mismatch !" << std::endl;
            exit(1);
        }
    }

    unsigned int jmin,jmax;
    if (col>0) {
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
    unsigned count = 0;

    if ((norm1>1e-4)&(norm2>1e-4)&(mat1.nlin()!=1)) {
        for(unsigned int i=0;i<mat1.nlin();++i) {
            for(unsigned int j=jmin;j<jmax;++j) {
                diff = std::abs(mat1(i,j)/norm1 - mat2(i,j)/norm2);
                flag = flag && (diff<eps);
                if (!(diff<eps)&&(++count<100)) {
                    if (count == 0)
                        std::cout << "ERROR NORM  mat(i,j) | mat2(i,j) | diff" << std::endl;
                    std::cout.precision(20);
                    std::cout << "ERROR NORM  " << mat1(i,j) << "  " << mat2(i,j) << "  " << diff << std::endl;
                    std::cout.flush();
                }
                if ( count >= 100 ) {
                    std::cout << "values from 1 to 100.... stopping display..." << std::endl;
                    break;
                }
            }
            if ( count >= 100 ) {
                break;
            }
        }
    } else {
        for(unsigned int i=0;i<mat1.nlin();++i) {
            for(unsigned int j=jmin;j<jmax;++j) {
                if (std::abs(mat2(i,j))>1e-4) {
                    diff = std::abs(mat1(i,j) - mat2(i,j))/std::abs(mat2(i,j));
                    flag = flag && (diff<eps);
                    if (!(diff<eps)&&(++count<100)) {
                        std::cout.precision(20);
                        std::cout << "ERROR RELATIVE  " << mat1(i,j) << "  " << mat2(i,j) << "  " << diff << std::endl;
                        std::cout.flush();
                    }
                    if ( count >= 100 ) {
                        std::cout << "values from 1 to 100.... stopping display..." << std::endl;
                        break;
                    }
                }
                else {
                    diff = std::abs(mat1(i,j) - mat2(i,j));
                    flag = flag && (diff<eps);
                    if (!(diff<eps)&&(++count<100)) {
                        std::cout.precision(20);
                        std::cout << "ERROR DIFF  " << mat1(i,j) << "  " << mat2(i,j) << "  " << diff << std::endl;
                        std::cout.flush();
                    }
                    if ( count >= 100 ) {
                        std::cout << "values from 1 to 100. Stopping the display..." << std::endl;
                        break;
                    }
                }
            }
            if ( count >= 100 ) {
                break;
            }
        }
    }
    return flag;
}

bool compare(maths::ifstream& ifs1,Vector& V1,maths::ifstream& ifs2,Vector& V2,const double eps) {

    try {
        ifs1 >> V1;
        ifs2 >> V2;
    } catch (std::string s) {
        std::cerr << s << std::endl;
        exit(1);
    }

    if (V1.size()!=V2.size()) {
        std::cerr << "ERROR : Dimension mismatch !" << std::endl;
        exit(1);
    }

    bool flag = true;

    const double norm1 = normInf(V1);
    const double norm2 = normInf(V2);

    if ((norm1>1e-4) && (norm2>1e-4)) {
        for (unsigned i=0;i<V1.size();++i) {
            const double diff = std::abs(V1(i)/norm1-V2(i)/norm2);
            flag = flag && (diff<eps);
            if (!(diff<eps)) {
                std::cout << "ERROR NORM  " << V1(i) << "  " << V2(i) << "  " << diff << std::endl;
                std::cout.flush();
            }
        }
    } else {
        for (unsigned i=0;i<V1.size();++i) {
            if (std::abs(V2(i))>1e-4) {
                const double diff = std::abs(V1(i)-V2(i))/std::abs(V2(i));
                flag = flag && (diff<eps);
                if (!(diff<eps)) {
                    std::cout << "ERROR RELATIVE  " << V1(i) << "  " << V2(i) << "  " << diff << std::endl;
                    std::cout.flush();
                }
            } else {
                const double diff = std::abs(V1(i)-V2(i));
                flag = flag && (diff<eps);
                if (!(diff<eps)) {
                    std::cout << "ERROR DIFF  " << V1(i) << "  " << V2(i) << "  " << diff << std::endl;
                    std::cout.flush();
                }
            }
        }
    }
    return flag;
}

template<class T>
bool compare_rdm(const T& mat1, const T& mat2, double eps, size_t col){
    // T is a Matrix

    if (col) {
        if ((mat1.ncol()<col) || (mat2.ncol()<col)) {
            std::cerr << "ERROR : Bad Column Id for matrices dimensions !" << std::endl;
            exit(1);
        }
    } else {
        if ((mat1.ncol()!=mat2.ncol()) || (mat1.nlin()!=mat2.nlin())) {
            std::cerr << "ERROR : Dimension mismatch !" << std::endl;
            exit(1);
        }
    }

    bool flag = true;
    double diff;
    unsigned int jmin,jmax;
    if (col>0) {
        jmin = col-1;
        jmax = col;
    } else {
        jmin = 0;
        jmax = mat1.ncol();;
    }

    std::cout.precision(20); // TODO
    for(unsigned int j=jmin;j<jmax;++j) {
        Vector col1 = mat1.getcol(j);
        Vector col2 = mat2.getcol(j);
        col1 = col1 - col1.mean();
        col2 = col2 - col2.mean();
        col1 = col1 / col1.norm();
        col2 = col2 / col2.norm();
        diff = (col1 - col2).norm();

        flag = flag && (diff<eps);
        if (diff>eps) {
            std::cout << "ERROR RDM (column " << j << " ) " << diff << std::endl;
            std::cout.flush();
        }
    }
    return flag;
}

template<class T>
bool compare_mag(const T& mat1, const T& mat2, double eps, size_t col){
    // T is a Matrix

    if (col) {
        if ((mat1.ncol()<col) || (mat2.ncol()<col)) {
            std::cerr << "ERROR : Bad Column Id for matrices dimensions !" << std::endl;
            exit(1);
        }
    } else {
        if ((mat1.ncol()!=mat2.ncol()) || (mat1.nlin()!=mat2.nlin())) {
            std::cerr << "ERROR : Dimension mismatch !" << std::endl;
            exit(1);
        }
    }

    bool flag = true;
    double diff;
    unsigned int jmin,jmax;
    if (col>0) {
        jmin = col-1;
        jmax = col;
    } else {
        jmin = 0;
        jmax = mat1.ncol();;
    }

    for(unsigned int j=jmin;j<jmax;++j) {
        Vector col1 = mat1.getcol(j);
        Vector col2 = mat2.getcol(j);
        col1 = col1 - col1.mean();
        col2 = col2 - col2.mean();
        diff = std::abs(1-col1.norm()/col2.norm()); //distance to 1

        flag = flag && (diff<eps);
        if (diff>eps) {
            std::cout << "ERROR MAG (column " << j << " ) = " << col1.norm()/col2.norm() << "\trelMAG = "<< diff << std::endl;
            std::cout.flush();
        }
    }
    return flag;
}
