// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

// Mesh simplification adapted from Kai Dang 2015 @ Inria

#include <mesh.h>
#include <options.h>
#include <cgal_lib.h>

using namespace OpenMEEG;

int main(int argc, char **argv)
{
    command_usage("Decimate a mesh:");
    const char * input_filename  = command_option("-i", (const char *) NULL,"Input image or mesh");
    const int    nb_points       = command_option("-n", 1000,"desired output number of vertices");
    const char * output_filename = command_option("-o", (const char *) NULL,"Output Mesh");

    Mesh m_in(input_filename, false);
    std::cout << "Input surface:\n nb of points: " << m_in.nb_vertices() << "\t nb of triangles:\t" << m_in.nb_triangles() << std::endl;
    Mesh m = cgal_decimate(m_in, nb_points);
    m.info();
    m.save(output_filename);

    return 0;
}
