// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <algorithm>

#include <constants.h>
#include <boundingbox.h>
#include <interface.h>

namespace OpenMEEG {

    /// Computes the total solid angle of a surface for a point p and tells whether p is inside the mesh or not.

    bool Interface::contains(const Vect3& p) const {
        const double solangle = solid_angle(p);

        if (almost_equal(solangle,-4*Pi))
            return true;

        if (almost_equal(solangle,0.0))
            return false;

        if (almost_equal(solangle,4*Pi)) {
            std::cerr << "Interface::contains_point(" << p << ") Error. This should not happen. Are you sure the mesh is properly oriented ?\n";
            return true;
        }

        std::cerr << "Interface::contains_point(" << p << ") Error. Are you sure the interface \"" << name() << "\" is closed? Solid angle: " << std::abs(solangle)/Pi <<"*PI." << std::endl;
        //  TODO: Is it really sane to return something ? Throw an exception instead...
        return (std::abs(solangle)>2*Pi) ? true : false;
    }

    /// Compute the solid angle which should be +/-4*Pi for a closed mesh if p is inside
    /// (the sign depends on the interface orientation), and 0 if p is outside.

    double Interface::solid_angle(const Vect3& p) const {
        double solangle = 0.0;
        for (const auto& omesh : oriented_meshes())
            solangle += omesh.orientation()*omesh.mesh().solid_angle(p);
        return solangle;
    }

    void Interface::set_to_outermost() {
        for (auto& omesh : oriented_meshes())
            omesh.mesh().outermost() = true;
        outermost_interface = true;
    }

    /// Check the global orientation: that the triangles are correctly oriented (outward-pointing normal)

    bool Interface::is_mesh_orientations_coherent(const bool doublecheck) {

        /// compute the bounding box:

        BoundingBox bb;
        for (const auto& omesh : oriented_meshes())
            for (const auto& vertex : omesh.mesh().vertices())
                bb.add(vertex);

        const Vect3& bbcenter = bb.center();

        double solangle = solid_angle(bbcenter);

        //  If the bounding box center is not inside the interface,
        //  try to test another point chosen randomly inside the bounding box, until
        //  its solid angle is significantly non zero.

        if (almost_equal(solangle,0.0))
            while (almost_equal(solangle,0.)) {
                const Vertex& V = bb.random_point();
                solangle = solid_angle(V);
            }

        if (almost_equal(solangle,4*Pi)) {
            // Reorient the interface.
            // TODO: With very little work OpenMEEG could be insensitive to global orientation of the interfaces.
            // and we could could remove this.

            std::cout << "Global reorientation of interface " << name() << std::endl;
            for (auto& omesh : oriented_meshes())
                omesh.change_orientation();
            solangle = -solangle;
        }

        if (almost_equal(solangle,0.0) || almost_equal(solangle,-4*Pi))
            return true;

        std::cout << solangle/Pi << "PI" << std::endl;

        //  In case of a bad random point location (too close to the mesh), do a double check:

        return doublecheck ? false : is_mesh_orientations_coherent(true);
    }
}
