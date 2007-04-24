#include "options.h"
#include "matrice.h"
#include "symmatrice.h"
#include "vecteur.h"
#include "om_utils.h"

int main( int argc, char** argv)
{
    command_usage("Convert squids in text file to a vtk file for vizualisation");
    const char *input_filename = command_option("-i",(const char *) SRCPATH("tools/data/MEGPositions.squids"),"Squids positions in original coordinate system");
    const char *output_filename = command_option("-o",(const char *) "MEGPositions.vtk","Squids positions with orientations in vtk format");
    if (command_option("-h",(const char *)0,0)) return 0;

    matrice squids(input_filename);
    assert(squids.nlin() == 151);

    FILE* f = fopen(output_filename,"w");
    if (f==NULL)
    {
        perror("fopen");
        return -1;
    }
    fprintf(f,"# vtk DataFile Version 3.0\n");
    fprintf(f,"vtk output\n");
    fprintf(f,"ASCII\n");
    fprintf(f,"DATASET POLYDATA\n");
    fprintf(f,"POINTS %d float\n",squids.nlin());
    for( unsigned int i = 0; i < squids.nlin(); i += 1 )
    {
        fprintf(f, "%f %f %f\n", squids(i,0), squids(i,1), squids(i,2));
    }
    fprintf(f,"POINT_DATA %d\n",squids.nlin());
    fprintf(f,"NORMALS normals float\n");
    for( unsigned int i = 0; i < squids.nlin(); i += 1 )
    {
        fprintf(f, "%f %f %f\n", squids(i,3), squids(i,4), squids(i,5));
    }
    fclose(f);
    return 0;
}
