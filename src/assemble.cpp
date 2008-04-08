/* FILE: $Id$ */

/*
Project Name : OpenMEEG

author            : $Author$
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
#include "fcontainer.h"
#include "cpuChrono.h"
#include "assemble.h"
#include "sensors.h"
#include "geometry.h"

using namespace std;

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
    * Computation of LHS for BEM Symmetric formulation
    **********************************************************************************************/
    if (!strcmp(argv[1],"-LHS")) {
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

        // Assembling matrix from discretization :
        LHS_matrice lhs(geo,GaussOrder);
        lhs.SAVE(argv[4]);
    }

    /*********************************************************************************************
    * Computation of general RHS from BEM Symmetric formulation
    **********************************************************************************************/
    else if (!strcmp(argv[1],"-RHS")) {

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
        bool checkClosedSurface = false;
        mesh_sources.load(argv[4],checkClosedSurface); // Load mesh without crashing when the surface is not closed

        // Assembling matrix from discretization :
        RHS_matrice mat(geo,mesh_sources,GaussOrder);
        mat.SAVE(argv[5]); // if outfile is specified
    }

    /*********************************************************************************************
    * Computation of RHS for discrete dipolar case
    **********************************************************************************************/
    else if(!strcmp(argv[1],"-rhsPOINT")) {

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

        // Loading matrix of dipoles :
        matrice dipoles(argv[4]);
        if(dipoles.ncol()!=6)
        {
            cerr << "Dipoles File Format Error" << endl;
            exit(1);
        }

        // Assembling matrix from discretization :
        unsigned int nd = (unsigned int) dipoles.nlin();
        std::vector<Vect3> Rs,Qs;
        for( unsigned int i=0; i<nd; i++ )
        {
            Vect3 r(3),q(3);
            for(int j=0;j<3;j++) r(j)   = dipoles(i,j);
            for(int j=3;j<6;j++) q(j-3) = dipoles(i,j);
            Rs.push_back(r); Qs.push_back(q);
        }
        RHSdip_matrice mat(geo, Rs, Qs, GaussOrder);

        // Saving RHS matrix for dipolar case :
        mat.SAVE(argv[5]);
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
	
	matrice source(newtaille,sourcetaille);
        matrice airescalp(newtaille,sourcetaille);
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
	int taille=geo.size();
        int sourcetaille = (geo.getM(geo.nb()-1)).nbTrgs();
	int newtaille=taille-sourcetaille;
	matrice source;
	source.loadBin(argv[4]);
	sparse_matrice stimelec;
        stimelec.loadBin(argv[5]);
	std::cout << "EITsource nlines: " << source.nlin() << "   ncols: " << source.ncol() << std::endl;
	std::cout << "stimelec nlines: " << stimelec.nlin() << "   ncols: " << stimelec.ncol() << std::endl;
        matrice stim(source.nlin(),stimelec.ncol());
	//        stim = source*stimelec;
	stim.saveBin(argv[6]);
    }


    /*********************************************************************************************
    * RK: Computation of RHS for discrete dipolar case: gradient wrt dipoles position and intensity!
    **********************************************************************************************/
    else if(!strcmp(argv[1],"-rhsPOINTgrad")) {

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

        // Loading matrix of dipoles :
        matrice dipoles(argv[4]);
        if(dipoles.ncol()!=6)
        {
            cerr << "Dipoles File Format Error" << endl;
            exit(1);
        }

        // Assembling matrix from discretization :
        unsigned int nd = (unsigned int) dipoles.nlin();
        std::vector<Vect3> Rs,Qs;
        for( unsigned int i=0; i<nd; i++ )
        {
            Vect3 r(3),q(3);
            for(int j=0;j<3;j++) r(j)   = dipoles(i,j);
            for(int j=3;j<6;j++) q(j-3) = dipoles(i,j);
            Rs.push_back(r); Qs.push_back(q);
        }

        RHSdip_grad_matrice mat( geo, Rs, Qs, GaussOrder);
        // Saving RHS matrix for dipolar case :
        mat.SAVE(argv[5]);
    }

    /*********************************************************************************************
    * Computation of the linear application which maps x (the unknown vector in symmetric system)
    * |----> v (potential at the electrodes)
    **********************************************************************************************/
    else if(!strcmp(argv[1],"-vToEEG")) {

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
        matrice patches(argv[4]);

        // Assembling matrix from discretization :
        // vToEEG is the linear application which maps x |----> v
        vToEEG_matrice mat(geo,patches);
        // Saving vToEEG matrix :
        mat.SAVE(argv[5]);
    }

    /*********************************************************************************************
    * Computation of the linear application which maps x (the unknown vector in symmetric system)
    * |----> bFerguson (contrib to MEG response)
    **********************************************************************************************/
    else if(!strcmp(argv[1],"-vToMEG")) {

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

        // Assembling matrix from discretization :
        vToMEG_matrice mat(geo,sensors);
        // Saving xToMEGrespCont matrix :
        mat.SAVE(argv[5]); // if outfile is specified
    }

    /*********************************************************************************************
    * Computation of the linear application which maps x (the unknown vector in symmetric system)
    * |----> binf (contrib to MEG response)
    **********************************************************************************************/
    else if(!strcmp(argv[1],"-sToMEG")) {

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
        bool checkClosedSurface = false;
        mesh_sources.load(argv[2],checkClosedSurface); // Load mesh without crashing when the surface is not closed

        // Load positions and orientations of sensors  :
        Sensors sensors(argv[3]);

        // Assembling matrix from discretization :
        sToMEG_matrice mat(mesh_sources, sensors);
        // Saving sToMEG matrix :
        mat.SAVE(argv[4]);
    }

    /*********************************************************************************************
    * Computation of the discrete linear application which maps x (the unknown vector in a symmetric system)
    * |----> binf (contrib to MEG response)
    **********************************************************************************************/
    // arguments are the positions and orientations of the squids,
    // the position and orientations of the sources and the output name.

    else if(!strcmp(argv[1],"-sToMEG_point")) {

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
        matrice dipoles(argv[2]);

        // Load positions and orientations of sensors  :
        Sensors sensors(argv[3]);

        sToMEGdip_matrice mat( dipoles, sensors );
        mat.SAVE(argv[4]);
    }
    /*********************************************************************************************
    * Computation of the discrete linear application which maps x (the unknown vector in a symmetric system)
    * |----> v, potential at a set of prescribed points within the 3D volume
    **********************************************************************************************/
    // arguments are the geom file
    // a tri-like file of point positions at which to evaluate the potential" << endl;

    else if(!strcmp(argv[1],"-SurfToVol")) {
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
	matrice points(argv[4]);
        SurfToVol_matrice mat(geo,points);
	std::cout << " mat(0,0) = " << mat(0,0) << std::endl;
	std::cout << " mat(1770,3445) = " << mat(1770,3445) << std::endl;
        // Saving SurfToVol matrix :
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
    cout << "   -LHS :   Compute LHS from BEM symmetric formulation." << endl;
    cout << "            Arguments :" << endl;
    cout << "               geometry file (.geom)" << endl;
    cout << "               conductivity file (.cond)" << endl;
    cout << "               output LHS matrix" << endl << endl;

    cout << "   -RHS :   Compute RHS from BEM symmetric formulation. " << endl;
    cout << "            Arguments :" << endl;
    cout << "               geometry file (.geom)" << endl;
    cout << "               conductivity file (.cond)" << endl;
    cout << "               mesh of sources (.tri .vtk .mesh .bnd)" << endl;
    cout << "               output RHS matrix" << endl << endl;

    cout << "   -rhsPOINT :   Compute RHS for discrete dipolar case. " << endl;
    cout << "            Arguments :" << endl;
    cout << "               geometry file (.geom)" << endl;
    cout << "               conductivity file (.cond)" << endl;
    cout << "               dipoles positions and orientations" << endl;
    cout << "               output RHS matrix" << endl << endl;

    cout << "   -EITsource :  Compute RHS for scalp current injection. " << endl;
    cout << "            Arguments :" << endl;
    cout << "               geometry file (.geom)" << endl;
    cout << "               conductivity file (.cond)" << endl;
    cout << "               output EITsource" << endl;
    cout << "               output airescalp" << endl << endl;

    cout << "   -EITstim :  Compute matrix directly mapping injected current values to EIT RHS. " << endl;
    cout << "            Arguments :" << endl;
    cout << "               geometry file (.geom)" << endl;
    cout << "               conductivity file (.cond)" << endl;
    cout << "               input EITsource" << endl;
    //    cout << "               input airescalp" << endl; 
    cout << "               input stimelec" << endl;
    cout << "               output EITstim" << endl << endl;


    cout << "   -vToEEG :   Compute the linear application which maps the potential" << endl;
    cout << "            on the scalp to the EEG electrodes"  << endl;
    cout << "            Arguments :" << endl;
    cout << "               geometry file (.geom)" << endl;
    cout << "               conductivity file (.cond)" << endl;
    cout << "               file containing the positions of EEG patches (.patches)" << endl;
    cout << "               output vToEEG matrice" << endl << endl;

    cout << "   -vToMEG :   Compute the linear application which maps the potential" << endl;
    cout << "            on the scalp to the MEG sensors"  << endl;
    cout << "            Arguments :" << endl;
    cout << "               geometry file (.geom)" << endl;
    cout << "               conductivity file (.cond)" << endl;
    cout << "               file containing the positions and orientations of the MEG sensors (.squids)" << endl;
    cout << "               output xToMEG matrix" << endl << endl;

    cout << "   -sToMEG :   Compute the linear application which maps the current" << endl;
    cout << "            dipoles on the source mesh to the MEG sensors" << endl;
    cout << "            Arguments :" << endl;
    cout << "               mesh file for distributed sources (.tri .vtk .mesh .bnd)" << endl;
    cout << "               positions and orientations of the MEG sensors (.squids)" << endl;
    cout << "               output sToMEG matrix" << endl << endl;

    cout << "   -sToMEG_point :   Compute the linear application which maps the current" << endl;
    cout << "            dipoles to the MEG sensors" << endl;
    cout << "            Arguments :" << endl;
    cout << "               dipoles positions and orientations" << endl;
    cout << "               positions and orientations of the MEG sensors (.squids)" << endl;
    cout << "               name of the output sToMEG matrix" << endl << endl;

    cout << "   -SurfToVol :   Compute the linear application which maps the surface potential" << endl;
    cout << "            and normal current to the value of the potential at a set of points in the volume" << endl;
    cout << "            Arguments :" << endl;
    cout << "               geom file" << endl;
     cout << "              cond file" << endl;
    cout << "               a tri file of point positions at which to evaluate the potential" << endl;
    cout << "               name of the output SurfToVol matrix" << endl << endl;


    exit(0);
}

