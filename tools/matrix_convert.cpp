#include "symmatrice.h"
#include "matrice.h"

#include "options.h"

using namespace std;

int main( int argc, char **argv)
{
    command_usage("Convert matrices between different formats");
    const char *input_filename = command_option("-i",(const char *) NULL,"Input Matrice");
    const char *output_filename = command_option("-o",(const char *) NULL,"Output Matrice");
    if (command_option("-h",(const char *)0,0)) return 0;

    matrice M;
    M.load(input_filename);
    M.save(output_filename);

    return 0;
}
