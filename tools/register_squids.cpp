/* FILE: $Id: register_squids.cpp 208 2008-02-29 13:28:33Z gramfort $ */

/*
Project Name : OpenMEEG

author            : $Author: gramfort $
version           : $Revision: 208 $
last revision     : $Date: 2008-02-29 14:28:33 +0100 (Fri, 29 Feb 2008) $
modified by       : $LastChangedBy: gramfort $
last modified     : $LastChangedDate: 2008-02-29 14:28:33 +0100 (Fri, 29 Feb 2008) $

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre 
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
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

#include "options.h"
#include "matrice.h"
#include "symmatrice.h"
#include "vecteur.h"
#include "om_utils.h"
#include "sensors.h"
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
    const double scale = command_option("-scale",10.0,"Scaling (10 by default for CTF data)"); // CTF counts in cm whereas MRI is in mm
    if (command_option("-h",(const char *)0,0)) return 0;

    Sensors squids;
    squids.load(squids_filename);
    size_t nb_positions = squids.getNumberOfPositions();

    matrice fiducials; fiducials.loadTxt(fiducials_filename);
    
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
    
    for( size_t i = 0; i < nb_positions; i += 1 )
    {
        vecteur position = squids.getPosition(i);
        vecteur orientation = squids.getOrientation(i);

        position = position*scale;
        position = R.transpose()*(position - T); // R is orthognal : R^-1 == R'
        orientation = R.inverse()*orientation;

        squids.setPosition(i,position);
        squids.setOrientation(i,orientation);
    }
    
    squids.save(squids_output_filename);
    if(rotation_filename != "") R.saveTxt(rotation_filename);
    if(translation_filename != "") T.saveTxt(translation_filename);
}
