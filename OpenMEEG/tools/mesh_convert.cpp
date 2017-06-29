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

#include <om_common.h>
#include <geometry.h>
#include <matrix.h>
#include <options.h>

using namespace std;
using namespace OpenMEEG;

//  Should not be here...

double determinant3x3(const Matrix& mat) {

    om_assert(mat.nlin() == mat.ncol());
    om_assert(mat.nlin() == 3);

    return mat(0,0)*(mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1))+
           mat(0,2)*(mat(1,0)*mat(2,1)-mat(1,1)*mat(2,0))+
           mat(0,1)*(mat(1,2)*mat(2,0)-mat(1,0)*mat(2,2));
}

int main( int argc, char **argv) {
    print_version(argv[0]);

    command_usage("Convert mesh between different formats");
    const char *input_filename = command_option("-i", (const char *) NULL, "Input Mesh");
    const char *output_filename = command_option("-o", (const char *) NULL, "Output Mesh");
    const double tx = command_option("-tx", 0.0, "Translation along the x axis");
    const double ty = command_option("-ty", 0.0, "Translation along the y axis");
    const double tz = command_option("-tz", 0.0, "Translation along the z axis");
    const double sx = command_option("-sx", 1.0, "Scaling along the x axis");
    const double sy = command_option("-sy", 1.0, "Scaling along the y axis");
    const double sz = command_option("-sz", 1.0, "Scaling along the z axis");
    const char* transfmat = command_option("-mat", (const char *) NULL, "4x4 Transformation Matrix (Assumed format ASCII)");
    const char* invert = command_option("-invert", (const char *) NULL, "Invert triangles point order");
    if (command_option("-h", (const char *)0, 0)) return 0;

    if(!input_filename || !output_filename) {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    Mesh m(input_filename);

    for ( Mesh::vertex_iterator vit = m.vertex_begin(); vit != m.vertex_end(); ++vit) {
        Vertex& v = **vit;
        v(0) = v(0)*sx + tx;
        v(1) = v(1)*sy + ty;
        v(2) = v(2)*sz + tz;
    }

    if ( transfmat ) {
        Matrix mat;
        mat.load(transfmat);

        if ((mat.nlin()!=4) || (mat.ncol()!=4)) {
            std::cerr << "Transformation matrix are 4x4 matrices, the one provided is " << mat.nlin() << 'x' << mat.ncol() << std::endl;
            return 2;
        }

        double mdet = determinant3x3(mat.submat(0, 3, 0, 3));
        if (mdet < 0 && !invert) // transformation is indirect => should force face flipping
            std::cout << "Warning : Transformation is indirect use -invert option to force face flipping" << std::endl;

        for ( Mesh::vertex_iterator vit = m.vertex_begin(); vit != m.vertex_end(); ++vit) {
            Vertex& v = **vit;
            Vector point(4);
            point.set(1.0);
            point(0) = v(0); point(1) = v(1); point(2) = v(2);
            Vector out_point = mat*point;
            v(0) = out_point(0); v(1) = out_point(1); v(2) = out_point(2);
        }
    }

    if ( invert ) {
        std::cout << "Running face flipping" << std::endl;
        m.flip_triangles();
    }

    m.correct_local_orientation();

    m.save(output_filename);

    return 0;
}
