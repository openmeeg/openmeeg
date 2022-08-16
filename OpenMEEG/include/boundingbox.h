// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

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
