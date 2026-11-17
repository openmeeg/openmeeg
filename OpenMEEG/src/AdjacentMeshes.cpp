// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <domain.h>
#include <AdjacentMeshes.h>
#include <OMExceptions.H>

namespace OpenMEEG {

    int AdjacentMeshes::compute_orientation() const {
        if (doms.empty())
            throw InvalidAdjacentMeshes(meshes[0]->name(),meshes[1]->name());
        if (meshes[0]==meshes[1])
            throw InvalidAdjacentMeshes(meshes[0]->name());

        // Adjacent meshes are meshes that share at least one domain with non null conductivity.
        // Mesh orientation should be coherent across domains so use only the first domain.

        return (doms[0]->mesh_orientation(*meshes[0])==doms[0]->mesh_orientation(*meshes[1])) ? 1 : -1;
    }
}
