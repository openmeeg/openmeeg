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

#include <interface.h>
#include <algorithm>
#include <ciso646>

namespace OpenMEEG {

    /// Computes the total solid angle of a surface for a point p and tells whether p is inside the mesh or not.
    bool Interface::contains_point(const Vect3& p) const 
    {
        double solangle = compute_solid_angle(p);

        if ( almost_equal(solangle, 0.) ) {
            return false;
        } else if ( almost_equal(solangle, -4.*M_PI) ) {
            return true;
        } else if ( almost_equal(solangle, 4.*M_PI) ) {
            std::cerr << "Interface::contains_point(" << p << ") Error. This should not happen. Are you sure the mesh is properly oriented ?\n";
            return true;
        } else {
            std::cerr << "Interface::contains_point(" << p << ") Error. Are you sure the interface \"" << name_ << "\" is closed? Solid angle: " << std::abs(solangle)/M_PI <<"*PI." << std::endl;
            return std::abs(solangle)>2*M_PI?true:false;
        }
    }

    /// compute the solid-angle which should be +/-4 * Pi for a closed mesh if p is inside, 0 if p is outside
    double Interface::compute_solid_angle(const Vect3& p) const 
    {
        double solangle = 0.0;
        for ( Interface::const_iterator omit = begin(); omit != end(); ++omit) {
                solangle += omit->orientation() * omit->mesh().compute_solid_angle(p);
        }
        return solangle;
    }

    void Interface::set_to_outermost() 
    {
        for ( Interface::iterator omit = begin(); omit != end(); ++omit) {
            omit->mesh().outermost() = true;
        }
        outermost_ = true;
    }

    /// Check the global orientation: that the triangles are correctly oriented (outward-pointing normal)
    bool Interface::check(bool checked)
    {
        /// compute the bounding box:
        double xmin = std::numeric_limits<double>::max();
        double ymin = std::numeric_limits<double>::max();
        double zmin = std::numeric_limits<double>::max();
        double xmax = -std::numeric_limits<double>::max();
        double ymax = -std::numeric_limits<double>::max();
        double zmax = -std::numeric_limits<double>::max();

        for ( Interface::const_iterator omit = begin(); omit != end(); ++omit) {
            for ( Mesh::const_vertex_iterator vit = omit->mesh().vertex_begin(); vit != omit->mesh().vertex_end(); ++vit) {
                xmin = std::min(xmin, (**vit).x());
                ymin = std::min(ymin, (**vit).y());
                zmin = std::min(zmin, (**vit).z());
                xmax = std::max(xmax, (**vit).x());
                ymax = std::max(ymax, (**vit).y());
                zmax = std::max(zmax, (**vit).z());
            }
        }
        
        Vect3 bbmin(xmin, ymin, zmin);
        Vect3 bbmax(xmax, ymax, zmax);
        Vect3 bbcenter = 0.5 * (bbmin + bbmax);

        // TODO if solidangle returns 0, then the center of BB is not inside, and thus we should randomly choose another inside point until solid_anlge gives + or -4 PI ?
        // else it causes a strange phenomenon for symmetric model, the point chosen for interface Cortex {north and south}, is on the surface...
        // compute the solid-angle from an inside point:
        double solangle = compute_solid_angle(bbcenter);
        bool closed;

        //if the bounding box center is not inside the interface,
        //we try to test another point inside the bounding box.
        if ( almost_equal(solangle, 0.)) {
            //std::cout<<"bbcenter is not inside interface: "<<name_<<std::endl;
            if ( not checked)
                std::srand((unsigned int)std::time(NULL));
            else
                std::srand((unsigned int)(std::time(NULL)+3583)); //the program runs faster than the change of time value

            while ( almost_equal(solangle, 0.) ) {
                Vect3 pt_rd((double)rand()/RAND_MAX*(xmax-xmin)+xmin,
                            (double)rand()/RAND_MAX*(ymax-ymin)+ymin,
                            (double)rand()/RAND_MAX*(zmax-zmin)+zmin);
                //std::cout<<"\ttest random point("<<pt_rd<<")\n";
                solangle=compute_solid_angle(pt_rd);
            }
        }

        if ( almost_equal(solangle, 0.) ) {
            closed = true;
        } else if ( almost_equal(solangle, -4.*M_PI)) {
            closed = true;
        } else if ( almost_equal(solangle, 4.*M_PI) ) {
            // TODO we still have to reorient the interface, but the code should be able with very little work to
            // be insensitive to global orientation of the interfaces, until that day, we reorient:
            std::cout << "Global reorientation of interface " << name() << std::endl;
            for ( Interface::iterator omit = begin(); omit != end(); ++omit) {
                omit->second = !omit->second;
            }
            closed = true;
        } else {
            std::cout << solangle/M_PI << "PI" << std::endl;
            //in case of a bad random point location (too close to the mesh), do a double check:
            closed = checked?false:this->check(true);
        }

        return closed;
    }
}
