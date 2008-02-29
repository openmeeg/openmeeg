/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
lastrevision      : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
*/

#include "mesh3.h"
#include "geometry.h"
#include "options.h"
#include <string>

using namespace std;

int main( int argc, char **argv)
{
    command_usage("Check mesh intersections in geometry file");
    const char *geom_filename = command_option("-g",(const char *) "","Input .geom file");
    const char *mesh_filename = command_option("-m",(const char *) NULL,"Mesh file (ex: to test .geom with cortex mesh)");
    if (command_option("-h",(const char *)0,0)) return 0;

    Geometry g;
    g.read(geom_filename);
    if (g.selfCheck()) {
        cout << ".geom : OK" << endl;
    }
    if(mesh_filename)
    {
        Mesh m;
        m.load(mesh_filename);
        if(g.check(m))
        {
            cout << ".geom and mesh : OK" << endl;
        }
    }
    return 0;
}
