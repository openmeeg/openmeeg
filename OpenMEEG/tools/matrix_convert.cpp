/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

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

#include <symmatrix.h>
#include <matrix.h>
#include <sparse_matrix.h>
#include <fast_sparse_matrix.h>
#include <fstream>
#include <options.h>

using namespace std;
using namespace OpenMEEG;

template <typename MATRIX>
void conversion(maths::ifstream& ifs,const char* input_format,maths::ofstream& ofs,const char* output_format,const char* output_filename) {
    MATRIX M;
    if (input_format)
        ifs >> maths::format(input_format) >> M;
    else
        ifs >> M;

    M.info();
    if (output_format)
        ofs << maths::format(output_format) << M;
    else
        ofs << maths::format(output_filename,maths::format::FromSuffix) << M;
}

int main(int argc, char **argv) try {

    print_version(argv[0]);

    command_usage("Convert full/sparse/symmetric vectors/matrices between different formats");
    const char* input_filename = command_option("-i",(const char *) NULL,"Input matrix/vector");
    const char* output_filename = command_option("-o",(const char *) NULL,"Output matrix/vector");
    const char* input_format = command_option("-if",(const char *) NULL,"Input file format : ascii, binary, tex, matlab");
    const char* output_format = command_option("-of",(const char *) NULL,"Output file format : ascii, binary, tex, matlab");
    if ( command_option("-h",(const char *)0,0) ) { return 0; }

    if ( argc < 2 || !input_filename || !output_filename ) {
        cout << "Not enough arguments, try the -h option" << endl;
        return 1;
    }

    maths::ifstream ifs(input_filename);
    maths::ofstream ofs(output_filename);

    try {
        conversion<Vector>(ifs,input_format,ofs,output_format,output_filename);
        return 0;
    } catch (OpenMEEG::maths::BadStorageType&) {
        //  Ignore storage type problems as they will tried in sequence.
    } catch (OpenMEEG::maths::BadVector&) {
        //  Bad vector may mean that this is not a vector but a matrix (handled below).
    } catch (OpenMEEG::maths::BadContent&) {
        //  Bad content type may mean that this another type of matrix/vector.
    } catch (...) {
        throw;
    }

    try {
        conversion<Matrix>(ifs,input_format,ofs,output_format,output_filename);
        return 0;
    } catch (OpenMEEG::maths::BadStorageType&) {
        //  Ignore storage type problems as they will tried in sequence.
    } catch (OpenMEEG::maths::BadContent&) {
        //  Bad content type may mean that this another type of matrix/vector.
    } catch (...) {
        throw;
    }

    try {
        conversion<SymMatrix>(ifs,input_format,ofs,output_format,output_filename);
        return 0;
    } catch (OpenMEEG::maths::BadStorageType&) {
        //  Ignore storage type problems as they will tried in sequence.
    } catch (OpenMEEG::maths::BadContent&) {
        //  Bad content type may mean that this another type of matrix/vector.
    } catch (...) {
        throw;
    }

    try {
        conversion<SparseMatrix>(ifs,input_format,ofs,output_format,output_filename);
        return 0;
    } catch (OpenMEEG::maths::BadStorageType&) {
        //  Ignore storage type problems as they will tried in sequence.
    } catch (OpenMEEG::maths::BadContent&) {
        //  Bad content type may mean that this another type of matrix/vector.
    } catch (...) {
        throw;
    }

    throw OpenMEEG::maths::ImpossibleObjectIdentification(input_filename);

} catch(OpenMEEG::maths::Exception& e) {
    std::cerr << e.what() << std::endl;
    return e.code();
} catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    return -1;
}
