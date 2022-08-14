// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <symmatrix.h>
#include <matrix.h>
#include <sparse_matrix.h>
#include <fast_sparse_matrix.h>
#include <fstream>
#include <commandline.h>

using namespace OpenMEEG;

template <typename MATRIX>
void conversion(maths::ifstream& ifs,const std::string& input_format,maths::ofstream& ofs,
                const std::string& output_format,const std::string& output_filename)
{
    MATRIX M;
    if (input_format!="")
        ifs >> maths::format(input_format) >> M;
    else
        ifs >> M;

    M.info();
    if (output_format!="")
        ofs << maths::format(output_format) << M;
    else
        ofs << maths::format(output_filename,maths::format::FromSuffix) << M;
}

int
main(int argc,char* argv[]) try {

    print_version(argv[0]);

    const CommandLine cmd(argc,argv,"Convert full/sparse/symmetric vectors/matrices between different formats");
    const std::string& input_filename  = cmd.option("-i", std::string(),"Input matrix/vector");
    const std::string& output_filename = cmd.option("-o", std::string(),"Output matrix/vector");
    const std::string& input_format    = cmd.option("-if",std::string(),"Input file format : ascii, binary, tex, matlab");
    const std::string& output_format   = cmd.option("-of",std::string(),"Output file format : ascii, binary, tex, matlab");

    if (cmd.help_mode())
        return 0;

    if (argc<2 || input_filename=="" || output_filename=="") {
        std::cout << "Missing arguments, try the -h option" << std::endl;
        return 1;
    }

    maths::ifstream ifs(input_filename.c_str());
    maths::ofstream ofs(output_filename.c_str());

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
