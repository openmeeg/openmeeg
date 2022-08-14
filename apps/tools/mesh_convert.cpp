// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <om_common.h>
#include <geometry.h>
#include <matrix.h>
#include <commandline.h>

using namespace OpenMEEG;

//  Should not be here...

double determinant3x3(const Matrix& mat) {

    om_assert(mat.nlin() == mat.ncol());
    om_assert(mat.nlin() == 3);

    return mat(0,0)*(mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1))+
           mat(0,2)*(mat(1,0)*mat(2,1)-mat(1,1)*mat(2,0))+
           mat(0,1)*(mat(1,2)*mat(2,0)-mat(1,0)*mat(2,2));
}

int
main(int argc,char* argv[]) {
    print_version(argv[0]);

    const CommandLine cmd(argc,argv,"Convert mesh between different formats");
    const std::string& input_filename  = cmd.option("-i",     std::string(),"Input Mesh");
    const std::string& output_filename = cmd.option("-o",     std::string(),"Output Mesh");
    const double       tx              = cmd.option("-tx",    0.0,          "Translation along the x axis");
    const double       ty              = cmd.option("-ty",    0.0,          "Translation along the y axis");
    const double       tz              = cmd.option("-tz",    0.0,          "Translation along the z axis");
    const double       sx              = cmd.option("-sx",    1.0,          "Scaling along the x axis");
    const double       sy              = cmd.option("-sy",    1.0,          "Scaling along the y axis");
    const double       sz              = cmd.option("-sz",    1.0,          "Scaling along the z axis");
    const std::string& transfmat       = cmd.option("-mat",   std::string(),"4x4 Transformation Matrix (ASCII format)");
    const bool         invert          = cmd.option("-invert",false,        "Invert triangles point order");

    if (cmd.help_mode())
        return 0;

    if (input_filename=="" || output_filename=="") {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    Mesh m(input_filename);

    for (const auto& vertex : m.vertices()) {
        Vertex& v = *vertex;
        v(0) = v(0)*sx+tx;
        v(1) = v(1)*sy+ty;
        v(2) = v(2)*sz+tz;
    }

    if (transfmat!="") {
        Matrix mat;
        mat.load(transfmat);

        if ((mat.nlin()!=4) || (mat.ncol()!=4)) {
            std::cerr << "Transformation matrix are 4x4 matrices, the one provided is " << mat.nlin() << 'x' << mat.ncol() << std::endl;
            return 2;
        }

        const double mdet = determinant3x3(mat.submat(0,3,0,3));
        if (mdet<0 && !invert) // transformation is indirect => should force face flipping
            std::cout << "Warning : Transformation is indirect use -invert option to force face flipping" << std::endl;

        for (const auto& vertex : m.vertices()) {
            Vertex& v = *vertex;
            Vector point(4);
            point.set(1.0);
            point(0) = v(0); point(1) = v(1); point(2) = v(2);
            Vector out_point = mat*point;
            v(0) = out_point(0); v(1) = out_point(1); v(2) = out_point(2);
        }
    }

    if (invert) {
        std::cout << "Change mesh orientation" << std::endl;
        m.change_orientation();
    }

    m.correct_local_orientation();

    m.save(output_filename);

    return 0;
}
