#include <matrice.h>
#include <symmatrice.h>
#include <vecteur.h>
#include "options.h"

using namespace std;

int main( int argc, char **argv)
{
    command_usage("Simple tool to select a few time frames in a dataset.");
    const char *input_filename = command_option("-i",(const char *) SRCPATH("data/Computations/Head1/Head1.src"),"Dataset from which frames are extracted");
    const char *output_filename = command_option("-o",(const char *) "extracted_frames.txt","Extracted time frames");
    const size_t first_frame = command_option("-f",0,"Index of first frame");
    const size_t length = command_option("-l",1,"Nb of frames to extract");
    const double mult = command_option("-m",1.0,"Nb of frames to extract");
    if (command_option("-h",(const char *)0,0)) return 0;

    matrice in(input_filename,'t');
    assert((first_frame+length) <= in.ncol());
    matrice out(in.nlin(),length);
    for( unsigned int i = 0; i < length; i += 1 )
    {
        out.setcol(i,in.getcol(first_frame+i));
    }
    if (mult != 1.0) {
        for (size_t j=0;j<out.ncol();j++)
            for (size_t i=0;i<out.nlin();i++)
            {
                out(i,j) = out(i,j)*mult;
            }
    }
    out.saveTxt(output_filename);
    return 0;
}
