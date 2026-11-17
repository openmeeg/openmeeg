// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#pragma once

#include <vector>
#include "mesh.h"

namespace OpenMEEG {

    class Domain;

    // Ideally, this should be defined in domain.h, but this would create a circular dependency,
    // so declare it here instead.

    using DomainsRefs = std::vector<const Domain*>;

    // Describes a pair of adjacent meshes, i.e. sharing at least one domain of non-null conductivity.
    // The two meshes must be different, otherwise an exception is triggered.

    class AdjacentMeshes {
    public:

        AdjacentMeshes(const Mesh& m1,const Mesh& m2,const DomainsRefs& d):
            meshes{&m1,&m2},doms(d),orientation(compute_orientation())
        { }

        const Mesh& operator()(const unsigned i) const { return *meshes[i]; }

        int relative_orientation() const { return orientation; }
        
        const DomainsRefs& common_domains() const { return doms; }

        // This is needed because MeshPair will be used to index a map.

        bool operator<(const AdjacentMeshes& mp1) const {
            if (meshes[0]==mp1.meshes[0])
                return meshes[1]<mp1.meshes[1];
            return meshes[0]<mp1.meshes[0];
        }

    private:

        int compute_orientation() const; // Cannot be inlined because Domains definition is not known here.

        const Mesh*       meshes[2];
        const DomainsRefs doms;
        const int         orientation;
    };

    using InteractingMeshes = std::vector<AdjacentMeshes>;
}
