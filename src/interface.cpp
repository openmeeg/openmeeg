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

    /**
     * Computes the total solid angle of a surface for a point p and tells whether p is inside the mesh or not.
     **/
    const bool Interface::contains_point(const Vect3& p) const 
    {
        double solangle = compute_solid_angle(p);

        if ( std::abs(solangle) < 1.e3*std::numeric_limits<double>::epsilon() ) {
            return false;
        } else if ( std::abs(solangle + 4.*M_PI) < 1.e3*std::numeric_limits<double>::epsilon()) {
            return true;
        } else if ( std::abs(solangle - 4.*M_PI) < 1.e3*std::numeric_limits<double>::epsilon()) {
            std::cerr << "Interface::contains_point(" << p << ") Error. Are you sure the mesh is properly oriented?\n";
            return false;
        } else {
            std::cerr << "Interface::contains_point(" << p << ") Error. Are you sure the mesh is closed?\n"
                << std::abs(solangle) << std::endl;
            exit(1);
        }
    }
    
    double Interface::compute_solid_angle(const Vect3& p) const // compute the solid-angle which should be +/- 4 * Pi for a closed mesh
    {
        double solangle = 0.0;
        for ( Interface::const_iterator omit = this->begin(); omit != this->end(); ++omit) {
            for ( Mesh::const_iterator tit = omit->mesh().begin(); tit != omit->mesh().end(); ++tit) {
                solangle += p.solangl((*tit).s1(), (*tit).s2(), (*tit).s3()) * omit->orientation();
            }
        }
        return solangle;
    }

    void Interface::set_to_outermost() 
    {
        for ( Interface::iterator omit = this->begin(); omit != this->end(); ++omit) {
            omit->mesh().outermost() = true;
        }
        outermost_ = true;
    }

    const bool Interface::closed() const
    {
        // compute the bounding box:
        double xmax = std::numeric_limits<double>::min();
        double ymax = std::numeric_limits<double>::min();
        double zmax = std::numeric_limits<double>::min();
        for ( Interface::const_iterator omit = this->begin(); omit != this->end(); ++omit) {
            for ( Mesh::const_vertex_iterator vit = omit->mesh().vertex_begin(); vit != omit->mesh().vertex_end(); ++vit) {
                xmax = std::max(xmax, (**vit).x());
                ymax = std::max(ymax, (**vit).y());
                zmax = std::max(zmax, (**vit).z());
            }
        }

        // compute the solid-angle from an outside point:
        double solangle = compute_solid_angle(2.*Vect3(xmax, ymax, zmax));

        if ( std::abs(solangle) < 1.e3*std::numeric_limits<double>::epsilon() ) {
            return true;
        }
        return false;
    }
}
