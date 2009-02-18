/* FILE: $Id$ */

/*
Project Name : OpenMEEG

version           : $Revision$
last revision     : $Date$
modified by       : $LastChangedBy$
last modified     : $LastChangedDate$

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
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

#define SAVEBIN

#ifdef SAVEBIN
    #define SAVE saveBin
#else
    #define SAVE saveTxt
#endif


#include <fstream>
#include <cstring>

#include "mesh3.h"
#include "integrator.h"
#include "cpuChrono.h"
#include "assemble.h"
#include "sensors.h"
#include "geometry.h"

using namespace std;
using namespace OpenMEEG;

int GaussOrder=3;
//int GaussOrder=0;

void getOutputFilepath(char* ref_filepath, char* output_filename, char* path);
void getHelp(char** argv);

int main(int argc, char** argv)
{
    if(argc < 2)
    {
        cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << endl;
        return 0;
    }

    if ((!strcmp(argv[1],"-h")) | (!strcmp(argv[1],"--help"))) getHelp(argv);

    disp_argv(argc,argv);

    // Start Chrono
    cpuChrono C;
    C.start();

    /*********************************************************************************************
    * Computation of Head Matrix for BEM Symmetric formulation
    **********************************************************************************************/
    if ((!strcmp(argv[1],"-HeadMat")) | (!strcmp(argv[1],"-HM")) | (!strcmp(argv[1],"-hm"))) {
        if (argc < 3)
        {
            std::cerr << "Please set geometry filepath !" << endl;
            exit(1);
        }
        if (argc < 4)
        {
            std::cerr << "Please set conductivities filepath !" << endl;
            exit(1);
        }
        if (argc < 5)
        {
            std::cerr << "Please set output filepath !" << endl;
            exit(1);
        }
        // Loading surfaces from geometry file
        Geometry geo;
        geo.read(argv[2],argv[3]);

        // Assembling Matrix from discretization :
        HeadMat HM(geo,GaussOrder);
        HM.SAVE(argv[4]);
    }

    /*********************************************************************************************
    * Computation of general Surface Source Matrix for BEM Symmetric formulation
    **********************************************************************************************/
    else if ((!strcmp(argv[1],"-SurfSourceMat")) | (!strcmp(argv[1],"-SSM")) | (!strcmp(argv[1],"-ssm"))) {
        if (argc < 3)
        {
            std::cerr << "Please set geometry filepath !" << endl;
            exit(1);
        }
        if (argc < 4)
        {
            std::cerr << "Please set conductivities filepath !" << endl;
            exit(1);
        }
        if (argc < 5)
        {
            std::cerr << "Please set 'mesh of sources' filepath !" << endl;
            exit(1);
        }

        // Loading surfaces from geometry file.
        Geometry geo;
        geo.read(argv[2],argv[3]);

        // Loading mesh for distributed sources
        Mesh  mesh_sources;
        mesh_sources.load(argv[4]);

        // Assembling Matrix from discretization :
        SurfSourceMat ssm(geo,mesh_sources,GaussOrder);
        ssm.SAVE(argv[5]); // if outfile is specified
    }

    /*********************************************************************************************
    * Computation of RHS for discrete dipolar case
    **********************************************************************************************/
    else if((!strcmp(argv[1],"-DipSourceMat")) |(!strcmp(argv[1],"-DSM"))|(!strcmp(argv[1],"-dsm")))  {
             if(argc < 3)
        {
            cerr << "Please set geometry filepath !" << endl;
            exit(1);
        }
        if (argc < 4)
        {
            std::cerr << "Please set conductivities filepath !" << endl;
            exit(1);
        }
        if(argc < 5)
        {
            cerr << "Please set dipoles filepath!" << endl;
            exit(1);
        }

        // Loading surfaces from geometry file.
        Geometry geo;
        geo.read(argv[2],argv[3]);

        // Loading Matrix of dipoles :
        Matrix dipoles(argv[4]);
        if(dipoles.ncol()!=6)
        {
            cerr << "Dipoles File Format Error" << endl;
            exit(1);
        }

        DipSourceMat dsm(geo, dipoles, GaussOrder);

        // Saving RHS Matrix for dipolar case :
        dsm.SAVE(argv[5]);
    }

    /*********************************************************************************************
    * Computation of RHS for EIT
    **********************************************************************************************/


    else if(!strcmp(argv[1],"-EITsource")) {

        if(argc < 3)
        {
            cerr << "Please set geometry filepath !" << endl;
            exit(1);
        }
        if (argc < 4)
        {
            std::cerr << "Please set conductivities filepath !" << endl;
            exit(1);
        }
        if (argc < 5)
        {
            std::cerr << "Please set output filepath !" << endl;
            exit(1);
        }
        if (argc < 6)
        {
            std::cerr << "Please set output filepath !" << endl;
            exit(1);
        }

        // Loading surfaces from geometry file.
        Geometry geo;
        geo.read(argv[2],argv[3]);

        int taille=geo.size();
        int sourcetaille = (geo.getM(geo.nb()-1)).nbTrgs();
        int newtaille=taille-sourcetaille;

        Matrix source(newtaille,sourcetaille);
        Matrix airescalp(newtaille,sourcetaille);
        source.set(0.0);
        airescalp.set(0.0);

        assemble_EITsource( geo, source, airescalp, GaussOrder);

        source.SAVE(argv[4]);
        airescalp.SAVE(argv[5]);
    }

    /*********************************************************************************************
    * Computation of RHS for EIT
    **********************************************************************************************/


    else if(!strcmp(argv[1],"-EITstim")) {

        if(argc < 3)
        {
            cerr << "Please set geometry filepath !" << endl;
            exit(1);
        }
        if (argc < 4)
        {
            std::cerr << "Please set conductivities filepath !" << endl;
            exit(1);
        }
        if (argc < 5)
        {
            std::cerr << "Please set EITsource filepath !" << endl;
            exit(1);
        }
        if (argc < 6)
        {
            std::cerr << "Please set stimelec filepath !" << endl;
            exit(1);
        }
        if (argc < 6)
        {
            std::cerr << "Please set output filepath !" << endl;
            exit(1);
        }
        // Loading surfaces from geometry file.
        Geometry geo;
        geo.read(argv[2],argv[3]);
        Matrix source;
        source.loadBin(argv[4]);
        SparseMatrix stimelec;
        stimelec.loadBin(argv[5]);
        Matrix stim(source.nlin(),stimelec.ncol());
        stim = source*stimelec;
        stim.saveBin(argv[6]);
    }

    /*********************************************************************************************
    * Computation of the linear application which maps the unknown vector in symmetric system,
    * (i.e. the potential and the normal current on all interfaces)
    * |----> v (potential at the electrodes)
    **********************************************************************************************/
    else if((!strcmp(argv[1],"-Head2EEGMat")) | (!strcmp(argv[1],"-H2EM"))| (!strcmp(argv[1],"-h2em"))) {

        if(argc < 3)
        {
            cerr << "Please set geometry filepath !" << endl;
            exit(1);
        }
        if (argc < 4)
        {
            std::cerr << "Please set conductivities filepath !" << endl;
            exit(1);
        }
        if(argc < 5)
        {
            cerr << "Please set patches filepath !" << endl;
            exit(1);
        }

        // Loading surfaces from geometry file.
        Geometry geo;
        geo.read(argv[2],argv[3]);

        // read the file containing the positions of the EEG patches
        Matrix patches(argv[4]);

        // Assembling Matrix from discretization :
        // Head2EEG is the linear application which maps x |----> v
        Head2EEGMat mat(geo,patches);
        mat.info();
        // Saving Head2EEG Matrix :
        mat.SAVE(argv[5]);
    }

    /*********************************************************************************************
    * Computation of the linear application which maps the unknown vector in symmetric system,
    * (i.e. the potential and the normal current on all interfaces)
    * |----> bFerguson (contrib to MEG response)
    **********************************************************************************************/
    else if((!strcmp(argv[1],"-Head2MEGMat"))|(!strcmp(argv[1],"-H2MM"))|(!strcmp(argv[1],"-h2mm"))) {

        if(argc < 3)
        {
            cerr << "Please set geometry filepath !" << endl;
            exit(1);
        }
        if (argc < 4)
        {
            std::cerr << "Please set conductivities filepath !" << endl;
            exit(1);
        }
        if(argc < 5)
        {
            cerr << "Please set squids filepath !" << endl;
            exit(1);
        }

        // Loading surfaces from geometry file.
        Geometry geo;
        geo.read(argv[2],argv[3]);

        // Load positions and orientations of sensors  :
        Sensors sensors(argv[4]);

        // Assembling Matrix from discretization :
        Head2MEGMat mat(geo,sensors);
        // Saving Head2MEG Matrix :
        mat.SAVE(argv[5]); // if outfile is specified
    }

    /*********************************************************************************************
    * Computation of the linear application which maps the distributed source
    * |----> binf (contrib to MEG response)
    **********************************************************************************************/
    else if((!strcmp(argv[1],"-SurfSource2MEGMat"))|(!strcmp(argv[1],"-SS2MM"))|(!strcmp(argv[1],"-ss2mm"))) {

        if(argc < 3)
        {
            cerr << "Please set 'mesh sources' filepath !" << endl;
            exit(1);
        }
        if(argc < 4)
        {
            cerr << "Please set squids filepath !" << endl;
            exit(1);
        }

        // Loading mesh for distributed sources :
        Mesh mesh_sources;
        mesh_sources.load(argv[2]);

        // Load positions and orientations of sensors  :
        Sensors sensors(argv[3]);

        // Assembling Matrix from discretization :
        SurfSource2MEGMat mat(mesh_sources, sensors);
        // Saving SurfSource2MEG Matrix :
        mat.SAVE(argv[4]);
    }

    /*********************************************************************************************
    * Computation of the discrete linear application which maps s (the dipolar source)
    * |----> binf (contrib to MEG response)
    **********************************************************************************************/
    // arguments are the positions and orientations of the squids,
    // the position and orientations of the sources and the output name.

    else if((!strcmp(argv[1],"-DipSource2MEGMat"))|(!strcmp(argv[1],"-DS2MM"))|(!strcmp(argv[1],"-ds2mm"))) {

        if (argc < 3)
        {
            cerr << "Please set dipoles filepath !" << endl;
            exit(1);
        }
        if (argc < 4)
        {
            cerr << "Please set squids filepath !" << endl;
            exit(1);
        }

        // Loading dipoles :
        Matrix dipoles(argv[2]);

        // Load positions and orientations of sensors  :
        Sensors sensors(argv[3]);

        DipSource2MEGMat mat( dipoles, sensors );
        mat.SAVE(argv[4]);
    }
    /*********************************************************************************************
    * Computation of the discrete linear application which maps x (the unknown vector in a symmetric system)
    * |----> v, potential at a set of prescribed points within the 3D volume
    **********************************************************************************************/
    // arguments are the geom file
    // a tri-like file of point positions at which to evaluate the potential" << endl;

    else if(!strcmp(argv[1],"-Surf2Vol")) {
        if (argc < 3)
        {
            cerr << "Please set geom filepath !" << endl;
            exit(1);
        }
        if (argc < 4)
        {
            cerr << "Please set cond filepath !" << endl;
            exit(1);
        }
        if (argc < 5)
        {
            cerr << "Please set point positions filepath !" << endl;
            exit(1);
        }
        if (argc < 6)
        {
            std::cerr << "Please set output filepath !" << endl;
            exit(1);
        }
        // Loading surfaces from geometry file
        Geometry geo;
        geo.read(argv[2],argv[3]);
        Matrix points(argv[4]);
        Surf2VolMat mat(geo,points);
        // Saving SurfToVol Matrix :
        mat.SAVE(argv[5]);
    }

    else cerr << "unknown argument: " << argv[1] << endl;

    // Stop Chrono
    C.stop();
    C.dispEllapsed();
}

void getOutputFilepath(char* ref_filepath, char* output_filename, char* path) {
    assert(path!=ref_filepath && path!=output_filename);
    // output filename on the same path as filename referenced in ref_filepath
    // go in search on all platform of path less filename included in ref_filepath.
    char* p = strrchr(ref_filepath, '/' );
    if (p == NULL)
        strcpy(path,output_filename);
    else
    {
        strncpy(path,ref_filepath,p-ref_filepath+1);
        strcat(path,output_filename);
    }
}

void getHelp(char** argv) {
    cout << argv[0] <<" [-option] [filepaths...]" << endl << endl;

    cout << "option :" << endl;
    cout << "   -HeadMat, -HM, -hm :   " << endl;
    cout << "    Compute Head Matrix for Symmetric BEM (left-hand side of linear system)." << endl;
    cout << "             Arguments :" << endl;
    cout << "               geometry file (.geom)" << endl;
    cout << "               conductivity file (.cond)" << endl;
    cout << "               output matrix" << endl << endl;

    cout << "   -SurfSourceMat, -SSM, -ssm :   " << endl;
    cout << "   Compute Surfacic Source Matrix for Symmetric BEM (right-hand side of linear system). " << endl;
    cout << "            Arguments :" << endl;
    cout << "               geometry file (.geom)" << endl;
    cout << "               conductivity file (.cond)" << endl;
    cout << "               mesh of sources (.tri .vtk .mesh .bnd)" << endl;
    cout << "               output matrix" << endl << endl;

    cout << "   -DipSourceMat, -DSM, -dsm:    " << endl;
    cout << "  Compute Dipolar Source Matrix for Symmetric BEM (right-hand side of linear system). " << endl;
    cout << "            Arguments :" << endl;
    cout << "               geometry file (.geom)" << endl;
    cout << "               conductivity file (.cond)" << endl;
    cout << "               dipoles positions and orientations" << endl;
    cout << "               output matrix" << endl << endl;

    cout << "   -EITsource :  Compute RHS for scalp current injection. " << endl;
    cout << "            Arguments :" << endl;
    cout << "               geometry file (.geom)" << endl;
    cout << "               conductivity file (.cond)" << endl;
    cout << "               output EITsource" << endl;
    cout << "               output matrix" << endl << endl;

    cout << "   -EITstim :  Compute Matrix directly mapping injected current values to EIT RHS. " << endl;
    cout << "            Arguments :" << endl;
    cout << "               geometry file (.geom)" << endl;
    cout << "               conductivity file (.cond)" << endl;
    cout << "               input EITsource" << endl;
    cout << "               input stimelec" << endl;
    cout << "               output matrix" << endl << endl;


    cout << "   -Head2EEGMat, -H2EM, -h2em : " << endl;
    cout << "        Compute the linear application which maps the potential" << endl;
    cout << "            on the scalp to the EEG electrodes"  << endl;
    cout << "            Arguments :" << endl;
    cout << "               geometry file (.geom)" << endl;
    cout << "               conductivity file (.cond)" << endl;
    cout << "               file containing the positions of EEG patches (.patches)" << endl;
    cout << "               output matrix" << endl << endl;

    cout << "   -Head2MEGMat, -H2MM, -h2mm : " << endl;
    cout << "          Compute the linear application which maps the potential" << endl;
    cout << "            on the scalp to the MEG sensors"  << endl;
    cout << "            Arguments :" << endl;
    cout << "               geometry file (.geom)" << endl;
    cout << "               conductivity file (.cond)" << endl;
    cout << "               file containing the positions and orientations of the MEG sensors (.squids)" << endl;
    cout << "               output matrix" << endl << endl;

    cout << "   -SurfSource2MEGMat, -SS2MM, -ss2mm : " << endl;
    cout << "         Compute the linear application which maps the " << endl;
    cout << "            distributed source  to the MEG sensors" << endl;
    cout << "            Arguments :" << endl;
    cout << "               mesh file for distributed sources (.tri .vtk .mesh .bnd)" << endl;
    cout << "               positions and orientations of the MEG sensors (.squids)" << endl;
    cout << "               output matrix" << endl << endl;

    cout << "   -DipSource2MEGMat, -DS2MM, -ds2mm :   Compute the linear application which maps the current" << endl;
    cout << "            dipoles to the MEG sensors" << endl;
    cout << "            Arguments :" << endl;
    cout << "               dipoles positions and orientations" << endl;
    cout << "               positions and orientations of the MEG sensors (.squids)" << endl;
    cout << "               name of the output matrix" << endl << endl;

    cout << "   -Surf2Vol :   Compute the linear application which maps the surface potential" << endl;
    cout << "            and normal current to the value of the potential at a set of points in the volume" << endl;
    cout << "            Arguments :" << endl;
    cout << "               geom file" << endl;
    cout << "               cond file" << endl;
    cout << "               a tri file of point positions at which to evaluate the potential" << endl;
    cout << "               name of the output matrix" << endl << endl;

    exit(0);
}
