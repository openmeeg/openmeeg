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

#include <interface.h>

namespace OpenMEEG {

    /// Computes the total solid angle of a surface for a point p and tells whether p is inside the mesh or not.
    const bool Interface::contains_point(const Vect3& p) const 
    {
        double solangle = compute_solid_angle(p);

        if ( std::abs(solangle) < 1.e3*std::numeric_limits<double>::epsilon() ) {
            return false;
        } else if ( std::abs(solangle + 4.*M_PI) < 1.e3*std::numeric_limits<double>::epsilon()) {
            return true;
        } else if ( std::abs(solangle - 4.*M_PI) < 1.e3*std::numeric_limits<double>::epsilon()) {
            std::cerr << "Interface::contains_point(" << p << ") Error. This should not happen. Are you sure the mesh is properly oriented ?\n";
            return false;
        } else {
            std::cerr << "Interface::contains_point(" << p << ") Error. Are you sure the interface \"" << name_ << "\" is closed?" << std::abs(solangle) << std::endl;
            exit(1);
        }
    }

    /// compute the solid-angle which should be +/-4 * Pi for a closed mesh if p is inside, 0 if p is outside
    double Interface::compute_solid_angle(const Vect3& p) const 
    {
        double solangle = 0.0;
        for ( Interface::const_iterator omit = begin(); omit != end(); ++omit) {
            for ( Mesh::const_iterator tit = omit->mesh().begin(); tit != omit->mesh().end(); ++tit) {
                solangle += p.solangl((*tit).s1(), (*tit).s2(), (*tit).s3())*omit->orientation();
            }
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
    const bool Interface::check()
    {
        /// compute the bounding box:
        double xmin = std::numeric_limits<double>::max();
        double ymin = std::numeric_limits<double>::max();
        double zmin = std::numeric_limits<double>::max();
        double xmax = std::numeric_limits<double>::min();
        double ymax = std::numeric_limits<double>::min();
        double zmax = std::numeric_limits<double>::min();

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

        // TODO if solidangle returns 0, then the center of BB is not inside, and thus we should randomly choose another inside point untill solid_anlge gives + or -4 PI ? else it causes a strange phenomenon for symmetric model, the point chosen for interface Cortex {north and south}, is on the surface...
        // compute the solid-angle from an inside point:
        double solangle = compute_solid_angle(bbcenter);
        bool closed;

        if ( std::abs(solangle) < 1.e3*std::numeric_limits<double>::epsilon() ) {
            closed = true;
        } else if ( std::abs(solangle + 4.*M_PI) < 1.e3*std::numeric_limits<double>::epsilon()) {
            closed = true;
        } else if ( std::abs(solangle - 4.*M_PI) < 1.e3*std::numeric_limits<double>::epsilon()) {
            std::cout << "Global Reorientation of interface " << name() << std::endl;
            for ( Interface::iterator omit = begin(); omit != end(); ++omit) {
                // omit->mesh().flip_triangles();
                omit->second = !omit->second; // TODO do we have to ?
            }
            closed = true;
        } else {
            closed = false;
        }

        return closed;
    }
}
