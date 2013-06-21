/*
Project Name : OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.sophia.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.sophia.inria.fr)

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

#include "geometry.h"
#include "options.h"

using namespace std;
using namespace OpenMEEG;

double determinant3x3(const Matrix& m);

int main( int argc, char **argv) {
    print_version(argv[0]);

    command_usage("Convert mesh between different formats");
    const char *input_filename = command_option("-i",(const char *) NULL,"Input Mesh");
    const char *output_filename = command_option("-o",(const char *) NULL,"Output Mesh");
    const double tx = command_option("-tx",0.0,"Translation along the x axis");
    const double ty = command_option("-ty",0.0,"Translation along the y axis");
    const double tz = command_option("-tz",0.0,"Translation along the z axis");
    const double sx = command_option("-sx",1.0,"Scaling along the x axis");
    const double sy = command_option("-sy",1.0,"Scaling along the y axis");
    const double sz = command_option("-sz",1.0,"Scaling along the z axis");
    const char* transfmat = command_option("-mat",(const char *) NULL,"4x4 Transformation Matrix (Assumed format ASCII)");
    const char* invert = command_option("-invert",(const char *) NULL,"Invert triangles point order");
    if (command_option("-h",(const char *)0,0)) return 0;

    if(!input_filename || !output_filename) {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    Mesh M;
    M.load_mesh(input_filename);

    for( int i = 0; i < M.nbPts(); ++i )
    {
        Vect3& pt = M[i];
        pt(0) = pt(0)*sx+tx;
        pt(1) = pt(1)*sy+ty;
        pt(2) = pt(2)*sz+tz;
    }

    if(transfmat)
    {
        Matrix m;
        m.load(transfmat);

        assert(m.nlin() == 4);
        assert(m.ncol() == 4);

        double mdet = determinant3x3(m.submat(0,3,0,3));
        if(mdet < 0 && !invert) // transformation is indirect => should force face flipping
        {
            cout << "Warning : Transformation is indirect use -invert option to force face flipping" << endl;
        }

        for( int i = 0; i < M.nbPts(); ++i )
        {
            Vect3& pt = M[i];
            Vector point(4);
            point.set(1.0);
            point(0) = pt(0); point(1) = pt(1); point(2) = pt(2);
            Vector out_point = m*point;
            pt(0) = out_point(0); pt(1) = out_point(1); pt(2) = out_point(2);

            Vect3& nm = M.normal(i);
            Vector normal(4);
            normal.set(0.0); // Hack to avoid the translation part
            normal(0) = nm(0); normal(1) = nm(1); normal(2) = nm(2);
            Vector out_normal = m*normal;
            nm(0) = out_normal(0); nm(1) = out_normal(1); nm(2) = out_normal(2);
        }
    }

    if(invert)
    {
        cout << "Running face flipping" << endl;
        M.flip_faces();
    }

    M.save(output_filename);

    return 0;
}

double determinant3x3(const Matrix& m) {
    assert(m.nlin() == m.ncol());
    assert(m.nlin() == 3);
    double f = 0.0;

    f += m(0,0)*m(1,1)*m(2,2);
    f += m(0,2)*m(1,0)*m(2,1);
    f += m(0,1)*m(1,2)*m(2,0);
    f -= m(0,2)*m(1,1)*m(2,0);
    f -= m(0,0)*m(1,2)*m(2,1);
    f -= m(0,1)*m(1,0)*m(2,2);

    return f;
}

