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

#pragma once

#include <random>

#include <vertex.h>

namespace OpenMEEG {

    /// An Oriented Mesh is a mesh associated with a boolean stating if it is well oriented.

    class BoundingBox {
    public:

        BoundingBox() { }

        void add(const Vertex& V) {
            xmin = std::min(xmin,V.x());
            ymin = std::min(ymin,V.y());
            zmin = std::min(zmin,V.z());
            xmax = std::max(xmax,V.x());
            ymax = std::max(ymax,V.y());
            zmax = std::max(zmax,V.z());
        }

        void add(const Vertex* Vp) { add(*Vp); }

        Vertex random_point() const {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> disx(xmin,xmax);
            std::uniform_real_distribution<> disy(ymin,ymax);
            std::uniform_real_distribution<> disz(zmin,zmax);
            return Vertex(disx(gen),disy(gen),disz(gen));
        }

        Vertex min() const { return Vertex(xmin,ymin,zmin); }
        Vertex max() const { return Vertex(xmax,ymax,zmax); }

        Vertex center() const { return 0.5*(min()+max()); }

    private:

        double xmin =  std::numeric_limits<double>::max();
        double ymin =  std::numeric_limits<double>::max();
        double zmin =  std::numeric_limits<double>::max();
        double xmax = -std::numeric_limits<double>::max();
        double ymax = -std::numeric_limits<double>::max();
        double zmax = -std::numeric_limits<double>::max();

    };
}
