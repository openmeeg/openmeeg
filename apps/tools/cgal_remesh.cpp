// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <mesh.h>
#include <options.h>
#include <cgal_lib.h>

using namespace OpenMEEG;


int main(int argc, char **argv) {
    command_usage("Re-mesh a mesh:");
    const char * input_filename  = command_option("-i",(const char *) NULL,"Input image or mesh");
    const double radius_bound    = command_option("-fs",1e-1,"facet radius bound of elements");
    const double distance_bound  = command_option("-fd",1e-1,"facet distance bound to the input surface");
    const char * output_filename = command_option("-o",(const char *) NULL,"Output Mesh");
    const char * sizing_field    = command_option("-field",(const char *) NULL,"(OPTIONAL) definition of the space to be refined 3 times finer (a matrix file: with either: \"x y z nx ny nz\" per line to define planes (by intersection of domains), or \"x y z r\" to define spheres (by union of domains).)");

    if ( command_option("-h",(const char *)0,0) ) {
        return 0;
    }
    if ( output_filename == NULL ) {
        std::cerr << "Set an output filename" << std::endl;
        return 0;
    }

    // Mesh input
    Mesh m_in(input_filename, false);
    std::cout << "Input surface:\n nb of points: " << m_in.nb_vertices() << "\t nb of triangles:\t" << m_in.nb_triangles() << std::endl;

    Mesh m_out = cgal_refine(m_in, radius_bound, distance_bound, sizing_field);
    m_out.save(output_filename);
    m_out.info();
    return 0;
}
