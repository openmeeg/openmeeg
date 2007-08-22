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
    command_usage("Convert squids positions from the machine coordinate system to the MRI coordinate system");
    const char *squids_filename = command_option("-i",(const char *) SRCPATH("tools/data/MEGPositions.squids"),"Squids positions in original coordinate system");
    const char *fiducials_orig_filename = command_option("-fo",(const char *) SRCPATH("tools/data/fiducials_orig.fid"),"Fiducial points in the original coordinate system (txt format or .hc CTF format)");
    const char *fiducials_dest_filename = command_option("-fd",(const char *) SRCPATH("tools/data/fiducials_dest.fid"),"Fiducial points in the destination coordinate system");
    const char *rotation_filename = command_option("-r",(const char *) "","Rotation matrix");
    const char *translation_filename = command_option("-t",(const char *) "","Translation vector");
    const char *squids_output_filename = command_option("-o",(const char *) "NewPositions.squids","Squids positions in the destination coordinate system");
    if (command_option("-h",(const char *)0,0)) return 0;

    matrice squids(squids_filename);

    matrice fiducials_orig(fiducials_orig_filename);
    matrice fiducials_dest(fiducials_dest_filename);

    size_t nb_squids = 151;
    assert(squids.nlin() == nb_squids);
    assert(squids.ncol() == 6);
    assert(fiducials_orig.nlin() == fiducials_orig.ncol());
    assert(fiducials_dest.nlin() == fiducials_dest.ncol());
    assert(fiducials_orig.nlin() == 3);
    assert(fiducials_dest.nlin() == 3);

    vecteur nas_orig = fiducials_orig.getlin(0);
    vecteur nas_dest = fiducials_dest.getlin(0);
    vecteur lpa_orig = fiducials_orig.getlin(1);
    vecteur lpa_dest = fiducials_dest.getlin(1);
    vecteur rpa_orig = fiducials_orig.getlin(2);
    vecteur rpa_dest = fiducials_dest.getlin(2);

    vecteur origin_orig = (lpa_orig+rpa_orig)/2.0;
    vecteur origin_dest = (lpa_dest+rpa_dest)/2.0;
    vecteur vx_orig = (nas_orig-origin_orig);
    vecteur vx_dest = (nas_dest-origin_dest);
    vecteur vz_orig = cross_product(vx_orig, lpa_orig-rpa_orig);
    vecteur vz_dest = cross_product(vx_dest, lpa_dest-rpa_dest);
    vecteur vy_orig = cross_product(vz_orig,vx_orig);
    vecteur vy_dest = cross_product(vz_dest,vx_dest);

    vx_orig = vx_orig/vx_orig.norm();
    vy_orig = vy_orig/vy_orig.norm();
    vz_orig = vz_orig/vz_orig.norm();

    vx_dest = vx_dest/vx_dest.norm();
    vy_dest = vy_dest/vy_dest.norm();
    vz_dest = vz_dest/vz_dest.norm();

    matrice rot_orig(3,3);
    matrice rot_dest(3,3);
    rot_orig.setcol(0,vx_orig);
    rot_orig.setcol(1,vy_orig);
    rot_orig.setcol(2,vz_orig);
    rot_dest.setcol(0,vx_dest);
    rot_dest.setcol(1,vy_dest);
    rot_dest.setcol(2,vz_dest);

    // Get the scaling factor since for now we have an isometry with R and T
    double scaling_factor = (lpa_dest - rpa_dest).norm() / (lpa_orig - rpa_orig).norm();
    scaling_factor += (lpa_dest - nas_dest).norm() / (lpa_orig - nas_orig).norm();
    scaling_factor += (rpa_dest - nas_dest).norm() / (rpa_orig - nas_orig).norm();
    scaling_factor = scaling_factor / 3.;

    matrice R =  rot_dest * rot_orig.inverse() * scaling_factor;
    matrice R_isometry =  rot_dest * rot_orig.inverse();
    vecteur T = origin_dest - R * origin_orig;

    for( unsigned int i = 0; i < nb_squids; i += 1 )
    {
        vecteur squid = squids.getlin(i);
        vecteur position(3);    vecteur orientation(3);
        position(0) = squid(0); orientation(0) = squid(3);
        position(1) = squid(1); orientation(1) = squid(4);
        position(2) = squid(2); orientation(2) = squid(5);
        position = R*position + T;
        orientation = R_isometry*orientation;
        squid(0) = position(0); squid(3) = orientation(0); // TODO : normalize direction
        squid(1) = position(1); squid(4) = orientation(1);
        squid(2) = position(2); squid(5) = orientation(2);
        squids.setlin(i,squid);
    }

    std::cout << "Storing new squids positions in : " << squids_output_filename << std::endl;
    squids.saveTxt(squids_output_filename);
    if(rotation_filename != "") R.saveTxt(rotation_filename);
    if(translation_filename != "") T.saveTxt(translation_filename);
}
