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

#include "interface.h"

namespace OpenMEEG {

    /**
     * Computes the total solid angle of a surface for a point p and tells whether p is inside the mesh or not.
     **/
    std::string Interface::keyword = "Mesh";

    bool Interface::contains_point(const Vect3& p) const {

        double solangle = 0.0;
        for (Interface::const_iterator mit = this->begin(); mit != this->end(); mit++) {
            for (Mesh::const_iterator tit = (*mit)->begin(); tit != (*mit)->end(); tit++) {
                solangle += p.solangl((*tit)(0).vertex(), (*tit)(1).vertex(), (*tit)(2).vertex());
            }
        }

        if (std::abs(solangle) < 1e3*std::numeric_limits<double>::epsilon()) {
            return false;
        } else if (std::abs(solangle + 4*M_PI) < 1e3*std::numeric_limits<double>::epsilon()) {
            return true;
        } else if (std::abs(solangle - 4*M_PI) < 1e3*std::numeric_limits<double>::epsilon()) {
            std::cerr << "Mesh::contains_point(" << p << ") Error. Are you sure the mesh is properly oriented?\n";
            return false;
        } else {
            std::cerr << "Mesh::contains_point(" << p << ") Error. Are you sure the mesh is closed?\n"
                << std::abs(solangle) << std::endl;
            exit(1);
        }
    }
    
    void  Interface::set_to_outermost() {
        for (PMeshes::iterator mit = meshes().begin(); mit != meshes().end(); mit++) {
            (*mit)->outermost() = true;
        }
        _outermost = true;
    }
}
