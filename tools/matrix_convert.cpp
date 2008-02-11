#include "symmatrice.h"
#include "matrice.h"
#include "sparse_matrice.h"
#include "fast_sparse_matrice.h"

#include "options.h"

using namespace std;

int main( int argc, char **argv) {
    command_usage("Convert matrices between different formats");
    const char *input_filename = command_option("-i",(const char *) NULL,"Input Matrice");
    const char *output_filename = command_option("-o",(const char *) NULL,"Output Matrice");
    const char *use_symmetric = command_option("-sym",(const char *) NULL,"Matrices are symmetric");
    const char *use_sparse = command_option("-sparse",(const char *) NULL,"Matrices are sparse");
    const char *use_binary = command_option("-bin",(const char *) NULL,"Input matrice is in binary format");
    const char *use_txt = command_option("-txt",(const char *) NULL,"Input matrices is in ascii format");
    const char *use_mat = command_option("-mat",(const char *) NULL,"Input matrices is in matlab format");
    if (command_option("-h",(const char *)0,0)) return 0;

    if (use_symmetric) {
        symmatrice M;
        if (use_binary) {
            M.loadBin(input_filename);
        } else if (use_txt) {
            M.loadTxt(input_filename);
        } else if (use_mat) {
            std::cerr << "Matlab format not supported for symmetric matrices" << std::endl;
        } else {
            M.load(input_filename);
        }
        M.save(output_filename);
    } else if (use_sparse) {
        sparse_matrice M;
        if (use_binary) {
            M.loadBin(input_filename);
        } else if (use_txt) {
            M.loadTxt(input_filename);
        } else if (use_mat) {
            M.loadMat(input_filename);
        } else {
            M.load(input_filename);
        }
        M.save(output_filename);
    } else {
        matrice M;
        if (use_binary) {
            M.loadBin(input_filename);
        } else if (use_txt) {
            M.loadTxt(input_filename);
        } else if (use_mat) {
            M.loadMat(input_filename);
        } else {
            M.load(input_filename);
        }
        M.save(output_filename);
    }

    return 0;
}
