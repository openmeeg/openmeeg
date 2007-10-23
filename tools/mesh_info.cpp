#include "mesh3.h"
#include "options.h"

using namespace std;

int main( int argc, char **argv)
{
    command_usage("Get info about a Mesh");
    const char *input_filename = command_option("-i",(const char *) "","Input Mesh");
    if (command_option("-h",(const char *)0,0)) return 0;

    Mesh* M = new Mesh();
    M->load(input_filename,false);
    M->save("toto.tri");

    delete M;
    return 0;
}
