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

#include "mesh.h"
#include "geometry.h"
#include "options.h"
#include <string>

using namespace OpenMEEG;

int main( int argc, char **argv)
{
    print_version(argv[0]);

    command_usage("Print Geometry information");
    const char *geom_filename = command_option("-g",(const char *) NULL,"Input .geom file");

    if (command_option("-h",(const char *)0,0)) return 0;

    if(!geom_filename) 
    {
        std::cout << "Not enough arguments, try the -h option" << std::endl;
        return 1;
    }

    int status = 0;
    Geometry geo;
    geo.read(geom_filename);

    if ( geo.selfCheck() )
    {
        std::cout << ".geom : OK" << std::endl;
    } else {
        status = 1;
    }
    geo.info();
    for ( Vertices::const_iterator vit = geo.vertex_begin(); vit != geo.vertex_end(); ++vit) {
        if ( (vit->index() == 69)||(vit->index() == 79)||(vit->index() == 95)||(vit->index() == 108)||(vit->index() == 112)||(vit->index() == 127)||(vit->index() == 600)||(vit->index() == 600)||(vit->index() == 600)||(vit->index() == 606)||(vit->index() == 606)||(vit->index() == 606)) {
            std::cout << "Index= " << vit->index() << "\t[" << *vit << "]" << std::endl;
        }
    }
    for ( Geometry::const_iterator mit = geo.begin(); mit != geo.end(); ++mit) {
        for ( Mesh::const_iterator tit = mit->begin(); tit != mit->end(); ++tit) {
            if ( (tit->index() == 69)||(tit->index() == 79)||(tit->index() == 95)||(tit->index() == 108)||(tit->index() == 112)||(tit->index() == 127)||(tit->index() == 600)||(tit->index() == 600)||(tit->index() == 600)||(tit->index() == 606)||(tit->index() == 606)||(tit->index() == 606)) {
                std::cout << "Triangle Index= " << tit->index() << "\t[" << tit->s1().index() << " , " << tit->s2().index() << " , " << tit->s3().index() << "]" << std::endl;
            }
        }
    }
    return status;
}
