#include "mesh3.h"
#include "options.h"

using namespace std;

int main( int argc, char **argv)
{
    command_usage("Convert mesh between different formats");
    const char *input_filename = command_option("-i",(const char *) NULL,"Input Mesh");
    const char *output_filename = command_option("-o",(const char *) NULL,"Output Mesh");
    const double tx = command_option("-tx",0.0,"Translation along the x axis");
    const double ty = command_option("-ty",0.0,"Translation along the y axis");
    const double tz = command_option("-tz",0.0,"Translation along the z axis");
    const double vx = command_option("-vx",1.0,"Scaling along the x axis");
    const double vy = command_option("-vy",1.0,"Scaling along the y axis");
    const double vz = command_option("-vz",1.0,"Scaling along the z axis");
    const bool apply_asa_flip = command_option("-flip",false,"Rotating axis if mesh comes from ASA");
    if (command_option("-h",(const char *)0,0)) return 0;

    Mesh M;
    M.load(input_filename,false);

    for( unsigned int i = 0; i < unsigned(M.nbPts()); i += 1 )
    {
        Vect3& pt = M[i];
        if (apply_asa_flip) {
            double tmp;
            tmp = pt(0);
            pt(0) = pt(1);
            pt(1) = tmp;
            pt(2) = -pt(2);
        }
        pt(0) = pt(0)+tx*vx;
        pt(1) = pt(1)+ty*vy;
        pt(2) = pt(2)+tz*vz;
    }

    M.save(output_filename);

    return 0;
}
