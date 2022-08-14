// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <MeshIOs/tri.h>
#include <MeshIOs/bnd.h>
#include <MeshIOs/off.h>
#include <MeshIOs/mesh.h>
#include <MeshIOs/vtk.h>
#include <MeshIOs/gifti.h>

namespace OpenMEEG {

    MeshIO::Registery MeshIO::registery;

    namespace MeshIOs {

        const Tri   Tri::prototype;
        const Bnd   Bnd::prototype;
        const Off   Off::prototype;
        const Mesh  Mesh::prototype;
        const Vtk   Vtk::prototype;
        const Gifti Gifti::prototype;
    }
}
