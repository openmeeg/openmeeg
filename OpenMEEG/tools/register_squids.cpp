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

#include "options.h"
#include "matrix.h"
#include "symmatrix.h"
#include "vector.h"
#include "om_utils.h"
#include "sensors.h"
#include <string>

using namespace OpenMEEG;

//  Should not be here.

Vector cross_product(const Vector &a, const Vector &b)
{
    om_assert(a.size() == 3);
    om_assert(b.size() == 3);
    Vector p(3);
    p(0) = a(1)*b(2)-a(2)*b(1);
    p(1) = a(2)*b(0)-a(0)*b(2);
    p(2) = a(0)*b(1)-a(1)*b(0);
    return p;
}

int main( int argc, char** argv)
{
    print_version(argv[0]);

    command_usage("Convert squids positions from the CTF MEG coordinate system to the MRI coordinate system");
    const char *squids_filename = command_option("-i",(const char *) NULL,"Squids positions in CTF coordinate system");
    const char *fiducials_filename = command_option("-f",(const char *) NULL,"Fiducial points in the MRI coordinate system (mm in txt format)");
    const char *rotation_filename = command_option("-r",(const char *) NULL,"Output Rotation Matrix");
    const char *translation_filename = command_option("-t",(const char *) NULL,"Output Translation vector");
    const char *squids_output_filename = command_option("-o",(const char *) NULL,"Squids positions in the MRI coordinate system");
    const double scale = command_option("-scale",10.0,"Scaling (10 by default for CTF data)"); // CTF counts in cm whereas MRI is in mm
    if (command_option("-h",(const char *)0,0)) return 0;

    if(!squids_filename || !fiducials_filename || !squids_output_filename) {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    Sensors squids;
    squids.load(squids_filename);
    size_t nb_positions = squids.getNumberOfPositions();

    Matrix fiducials; fiducials.load(fiducials_filename);

    if ((fiducials.nlin()!=3) ||
        (fiducials.ncol()!=3))
        throw std::invalid_argument("OpenMEEG only handles 3 3D fiducial points.");

    Vector nas = fiducials.getlin(0); // Nasion
    Vector lpa = fiducials.getlin(1); // Left preauricular
    Vector rpa = fiducials.getlin(2); // Right preauricular

    Vector origin = (lpa+rpa)/2.0;
    Vector vx = (nas-origin);
    Vector vz = cross_product(vx, lpa-rpa);
    Vector vy = cross_product(vz,vx);

    vx = vx/vx.norm();
    vy = vy/vy.norm();
    vz = vz/vz.norm();

    Matrix R(3,3);
    R.setlin(0,vx);
    R.setlin(1,vy);
    R.setlin(2,vz);

    Vector T = R * origin;
    T = T * (-1);

    for( size_t i = 0; i < nb_positions; i += 1 )
    {
        Vector position = squids.getPosition(i);
        Vector orientation = squids.getOrientation(i);

        position = position*scale;
        position = R.transpose()*(position - T); // R is orthognal : R^-1 == R'
        orientation = R.inverse()*orientation;

        squids.setPosition(i,position);
        squids.setOrientation(i,orientation);
    }

    squids.save(squids_output_filename);
    if(rotation_filename) R.save(rotation_filename);
    if(translation_filename) T.save(translation_filename);
}
