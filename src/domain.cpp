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

#include "domain.h"

namespace OpenMEEG {

    // Domain::Domain(const Domain& d) { *this = d; }

    // Domain& Domain::operator= (const Domain& d) {
        // if (this != &d) {
            // copy(d);
        // }
        // return *this;
    // }

    // void Domain::copy(const Domain& d) {
        // base::operator=(d);
        // sigma_     = d.sigma_;
        // outermost_ = d.outermost_;
        // innermost_ = d.innermost_;
    // }

    bool Domain::contains_point(const Vect3& p) const {
        bool inside = true;
        for (Domain::const_iterator hit = this->begin(); hit != this->end(); hit++) {
            inside = (inside &&  (hit->interface().contains_point(p) == hit->inside()));
        }
        return inside;
    }

    bool Domain::operator==(const Domain& d) const {
        // for (Domain::const_iterator dit1 = d1.begin(), dit2 = d2.begin(); dit1 != d1.end(); ++dit1, ++dit2) { TODO
            // if (&dit1->interface() != &dit2->interface()) {
                // return false;
            // }
        // }
        return true;
    }

    void Domain::destroy() {
        this->clear();
    }

    void Domain::info() const {

        std::cout << "Info:: Domain name : "  << name() << std::endl;
        std::cout << "\t\tConductivity : "    << sigma() << std::endl;
        std::cout << "\t\tComposed by interfaces : ";
        for (base::const_iterator hit = this->begin(); hit != this->end(); hit++) {
            if (hit->inside()) {
                std::cout << "-";
            } else {
                std::cout << "+";
            }
            std::cout << hit->interface().name() << " ";
        }
        std::cout << std::endl;
        if (innermost()) {
            std::cout << "\t\tConsidered as the innermost domain." << std::endl;
        }
        if (outermost()) {
            std::cout << "\t\tConsidered as the outermost domain." << std::endl;
        }
        for (base::const_iterator hit = this->begin(); hit != this->end(); hit++) {
            std::cout << "\t\tInterface \"" << hit->first.name() << "\"= { ";
            for (Interface::const_iterator mit = hit->first.begin(); mit != hit->first.end(); mit++) {
                std::cout << "Mesh \""<< (*mit)->name() << "\"";
                if ((*mit)->outermost()) {
                    std::cout << "(OUTERMOST)";
                }
                std::cout << ", ";
            }
            std::cout << " }" << std::endl;
        }
    }
}
