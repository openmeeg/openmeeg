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
    const char *output_filename = command_option("-o",(const char *) "","Output Mesh");
    const double smoothing_intensity = command_option("-s",0.1,"Smoothing Intensity");
    const size_t niter = command_option("-n",1000,"Number of iterations");
    if (command_option("-h",(const char *)0,0)) return 0;

    Mesh* M = new Mesh();
    M->load(input_filename);
    M->smooth(smoothing_intensity,niter);
    M->save(output_filename);

    delete M;
    return 0;
}
