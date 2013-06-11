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

#ifndef OPENMEEG_INTERFACE_H
#define OPENMEEG_INTERFACE_H

#include <string>
#include <vector>
#include "mesh.h"

namespace OpenMEEG {

    /** \brief  Interface

        Interface class

    **/

    class OPENMEEG_EXPORT Interface: protected std::vector<Mesh *> {

    public:

        // static void set_file_format(Mesh& mesh, const std::string& name) { 
            // mesh.get_file_format(name.c_str()); 
            // this->name = name;
        // }
        typedef std::vector<Mesh *> base;

        std::string keyword;

        Interface(std::string name): base(), keyword("Interface") {
            // this->push_back(&Mesh(name)); // TODO error: taking address of temporary
        }

        inline std::string   name() const   { return _name; }
        inline std::string & name()         { return _name; } // TODO maybe remove

        // inline std::istream& operator>>(std::istream &is, Mesh &m) {
            // if (strcmp(keyword, "Mesh") == 0) {
            // }
        // }

    private:
        std::string _name; // might be "i0" by default
    };

}

#endif  //  ! OPENMEEG_INTERFACE_H

