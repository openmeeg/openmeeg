/* FILE: $Id$ */

/*
Project Name : $Project$

author            : $Author$
version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

$License$
*/

#include "mesh3.h"
#include "options.h"

using namespace std;

int main( int argc, char **argv)
{
    command_usage("Concat 2 mesh and save the result");
    const char *input_filename1 = command_option("-i1",(const char *) "","Input Mesh 1");
    const char *input_filename2 = command_option("-i2",(const char *) "","Input Mesh 2");
    const char *output_filename = command_option("-o",(const char *) "","Output Mesh");
    if (command_option("-h",(const char *)0,0)) return 0;

    Mesh* M1 = new Mesh();
    M1->load(input_filename1);

    Mesh* M2 = new Mesh();
    M2->load(input_filename2);

    M1->append(M2);
    M1->save(output_filename);

    delete M1;
    delete M2;

    return 0;
}
