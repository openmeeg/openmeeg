// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <domain.h>

namespace OpenMEEG {

    bool Domain::contains(const Vect3& p) const {
        bool inside = true;
        for (const auto& boundary : boundaries())
            inside = (inside && (boundary.interface().contains(p)==boundary.inside()));
        return inside;
    }

    void Domain::info(const bool outermost) const {

        std::cout << "Info:: Domain name : "  << name() << std::endl;
        std::cout << "\t\tConductivity : "    << conductivity() << std::endl;
        std::cout << "\t\tComposed by interfaces : ";
        for (const auto& boundary : boundaries()) {
            std::cout << ((boundary.inside()) ? '-' : '+');
            std::cout << boundary.interface().name() << " ";
        }
        std::cout << std::endl;
        if (outermost)
            std::cout << "\t\tConsidered as the outermost domain." << std::endl;
        
        for (const auto& boundary : boundaries()) {
            std::cout << "\t\tInterface \"" << boundary.interface().name() << "\"= { ";
            for (const auto& oriented_mesh : boundary.interface().oriented_meshes()) {
                std::cout << "mesh \""<< oriented_mesh.mesh().name() << "\"";
                if (oriented_mesh.mesh().outermost())
                    std::cout << "(outermost)";
                std::cout << ", ";
            }
            std::cout << "\b\b }" << std::endl;
        }
    }
}
