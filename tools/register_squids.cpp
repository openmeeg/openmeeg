#include "options.h"
#include "matrice.h"
#include "symmatrice.h"
#include "vecteur.h"
#include "om_utils.h"
#include <string>

vecteur cross_product(const vecteur &a, const vecteur &b)
{
    assert(a.size() == 3);
    assert(b.size() == 3);
    vecteur p(3);
    p(0) = a(1)*b(2)-a(2)*b(1);
    p(1) = a(2)*b(0)-a(0)*b(2);
    p(2) = a(0)*b(1)-a(1)*b(0);
    return p;
}

int main( int argc, char** argv)
{
    command_usage("Convert squids positions from the CTF MEG coordinate system to the MRI coordinate system");
    const char *squids_filename = command_option("-i",(const char *) SRCPATH("tools/data/MEGPositions.squids"),"Squids positions in CTF coordinate system");
    const char *fiducials_filename = command_option("-f",(const char *) SRCPATH("tools/data/fiducials_orig.fid"),"Fiducial points in the MRI coordinate system (mm in txt format)");
    const char *rotation_filename = command_option("-r",(const char *) "","Output Rotation matrix");
    const char *translation_filename = command_option("-t",(const char *) "","Output Translation vector");
    const char *squids_output_filename = command_option("-o",(const char *) "NewPositions.squids","Squids positions in the MRI coordinate system");
    if (command_option("-h",(const char *)0,0)) return 0;

    matrice squids; squids.loadTxt(squids_filename);
    matrice fiducials; fiducials.loadTxt(fiducials_filename);

    size_t nb_squids = 151;
    assert(squids.nlin() == nb_squids);
    assert(squids.ncol() == 6);
    assert(fiducials.nlin() == 3);
    assert(fiducials.ncol() == 3);

    vecteur nas = fiducials.getlin(0); // Nasion
    vecteur lpa = fiducials.getlin(1); // Left preauricular
    vecteur rpa = fiducials.getlin(2); // Right preauricular

    vecteur origin = (lpa+rpa)/2.0;
    vecteur vx = (nas-origin);
    vecteur vz = cross_product(vx, lpa-rpa);
    vecteur vy = cross_product(vz,vx);

    vx = vx/vx.norm();
    vy = vy/vy.norm();
    vz = vz/vz.norm();

    matrice R(3,3);
    R.setlin(0,vx);
    R.setlin(1,vy);
    R.setlin(2,vz);

    vecteur T = R * origin;
    T = T * (-1);

    // std::cout << "Number of squids : " << nb_squids << std::endl;
    for( unsigned int i = 0; i < nb_squids; i += 1 )
    {
        vecteur squid = squids.getlin(i);
        vecteur position(3);    vecteur orientation(3);
        position(0) = squid(0); orientation(0) = squid(3);
        position(1) = squid(1); orientation(1) = squid(4);
        position(2) = squid(2); orientation(2) = squid(5);
        
        position = position*10; // CTF counts in cm whereas MRI is in mm
        position = R.transpose()*(position - T); // R is orthognal : R^-1 == R'
        orientation = R.inverse()*orientation;
        
        squid(0) = position(0); squid(3) = orientation(0);
        squid(1) = position(1); squid(4) = orientation(1);
        squid(2) = position(2); squid(5) = orientation(2);
        squids.setlin(i,squid);
    }

    // std::cout << "Storing new squids positions in : " << squids_output_filename << std::endl;
    squids.saveTxt(squids_output_filename);
    if(rotation_filename != "") R.saveTxt(rotation_filename);
    if(translation_filename != "") T.saveTxt(translation_filename);
}
