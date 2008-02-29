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
#include "options.h"

using namespace std;

int main( int argc, char **argv)
{
    command_usage("Get info about a Mesh");
    const char *input_filename = command_option("-i",(const char *) "","Input Mesh");
    if (command_option("-h",(const char *)0,0)) return 0;

    Mesh M;
    M.load(input_filename,false);

    return 0;
}
