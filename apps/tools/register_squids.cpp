// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <string>

#include "matrix.h"
#include "symmatrix.h"
#include "vector.h"
#include "sensors.h"

#include "commandline.h"

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

int
main(int argc,char* argv[]) {

    print_version(argv[0]);

    const CommandLine cmd(argc,argv,"Convert squids positions from the CTF MEG coordinate system to the MRI coordinate system");
    const std::string& squids_filename        = cmd.option("-i",std::string(),"Squids positions in CTF coordinate system");
    const std::string& fiducials_filename     = cmd.option("-f",std::string(),"Fiducial points in the MRI coordinate system (mm in txt format)");
    const std::string& rotation_filename      = cmd.option("-r",std::string(),"Output Rotation Matrix");
    const std::string& translation_filename   = cmd.option("-t",std::string(),"Output Translation vector");
    const std::string& squids_output_filename = cmd.option("-o",std::string(),"Squids positions in the MRI coordinate system");
    const double       scale                  = cmd.option("-scale",10.0,     "Scaling (10 by default for CTF data)"); // CTF uses cm whereas MRI is in mm

    if (cmd.help_mode())
        return 0;

    if(squids_filename=="" || fiducials_filename=="" || squids_output_filename=="") {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    Sensors squids;
    squids.load(squids_filename);
    const size_t nb_positions = squids.getNumberOfPositions();

    Matrix fiducials; fiducials.load(fiducials_filename);

    if ((fiducials.nlin()!=3) ||
        (fiducials.ncol()!=3))
        throw std::invalid_argument("OpenMEEG only handles 3 3D fiducial points.");

    const Vector nas = fiducials.getlin(0); // Nasion
    const Vector lpa = fiducials.getlin(1); // Left preauricular
    const Vector rpa = fiducials.getlin(2); // Right preauricular

    const Vector origin = (lpa+rpa)/2.0;
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

    const Vector& T = -(R*origin);

    for (size_t i=0; i<nb_positions; ++i) {

        Vector position = squids.getPosition(i);
        Vector orientation = squids.getOrientation(i);

        position = position*scale;
        position = R.transpose()*(position - T); // R is orthognal : R^-1 == R'
        orientation = R.inverse()*orientation;

        squids.setPosition(i,position);
        squids.setOrientation(i,orientation);
    }

    squids.save(squids_output_filename);

    if (rotation_filename!="")
        R.save(rotation_filename);

    if (translation_filename!="")
        T.save(translation_filename);

    return 0;
}
