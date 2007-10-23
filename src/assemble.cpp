#define SAVEBIN

#ifdef SAVEBIN
    #define SAVE saveBin
    #define SAVESUB saveSubBin
#else
    #define SAVE saveTxt
    #define SAVESUB saveSubTxt
#endif

#include <fstream>
#include <cstring>

#include "mesh3.h"
#include "integrator.h"
#include "fcontainer.h"
#include "cpuChrono.h"
#include "assemble.h"
#include "sensors.h"

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

    cout << endl << "| ------ " << argv[0] << " -------" << endl;
    for( int i = 1; i < argc; i += 1 )
    {
        cout << "| " << argv[i] << endl;
    }
    cout << "| -----------------------" << endl;

    // Start Chrono
    cpuChrono C;
    C.start();

    /*********************************************************************************************
    * Computation of LHS for BEM Symmetric formulation
    **********************************************************************************************/
    if (!strcmp(argv[1],"-LHS"))
    {
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
        // Loading surfaces from geometry file. 'taille' = sum on surfaces of number of points and number of triangles
        Geometry geo;
        unsigned int taille = geo.read(argv[2],argv[3]); // FIXME : read returns an int not unsigned

        // Assembling matrix from discretization :
        symmatrice mat(taille);
        assemble_LHS(geo,mat,GaussOrder);

        // Deflation the last diagonal bloc of new 'mat'  :
        int newtaille=taille-(geo.getM(geo.nb()-1)).nbTrgs();
        int offset=newtaille-(geo.getM(geo.nb()-1)).nbPts();
        deflat(mat,offset,newtaille-1,mat(offset,offset)/(newtaille-offset));

        mat.SAVESUB(argv[4],0,newtaille-1,0,newtaille-1);
    }

    /*********************************************************************************************
    * Computation of general RHS from BEM Symmetric formulation
    **********************************************************************************************/
    else if (!strcmp(argv[1],"-RHS"))
    {

        if (argc < 3)
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
            cerr << "Please set 'mesh of sources' filepath !" << endl;
            exit(1);
        }

        // Loading surfaces from geometry file.
        Geometry geo;
        int taille = geo.read(argv[2],argv[3]);

        // Loading mesh for distributed sources
        Mesh mesh_sources;
        bool checkClosedSurface = false;
        mesh_sources.load(argv[4],false); // Load mesh without crashing when the surface is not closed

        // Assembling matrix from discretization :
        int newtaille = taille-(geo.getM(geo.nb()-1)).nbTrgs();
        matrice mat(newtaille,mesh_sources.nbPts());
        mat.set(0.0);
        assemble_RHS(geo,mesh_sources,mat,GaussOrder);

        mat.SAVE(argv[5]); // if outfile is specified
    }

    /*********************************************************************************************
    * Computation of RHS from BEM Symmetric formulation with precomputation
    **********************************************************************************************/
    else if(!strcmp(argv[1],"-RHS2"))
    {

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
            cerr << "Please set 'mesh of sources' filepath !" << endl;
            exit(1);
        }

        // Loading surfaces from geometry file.
        Geometry geo;

        // Loading mesh for distributed sources
        Mesh mesh_sources;
        bool checkClosedSurface = false;
        mesh_sources.load(argv[4],false); // Load mesh without crashing when the surface is not closed

        // Assembling matrix from discretization :
        int newtaille = geo.getM(0).nbPts()+geo.getM(0).nbTrgs();
        matrice mat(newtaille,mesh_sources.nbPts()+mesh_sources.nbTrgs());
        mat.set(0.0);
        assemble_RHS2(geo,mesh_sources,mat,GaussOrder);

        // Saving RHS matrix :
        mat.SAVESUB(argv[5],0,newtaille-1,0,mesh_sources.nbPts()-1);
    }

    /*********************************************************************************************
    * Computation of RHS for discrete dipolar case
    **********************************************************************************************/
    else if(!strcmp(argv[1],"-rhsPOINT")){

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
        int taille=geo.read(argv[2],argv[3]);

        // Loading matrix of dipoles :
        matrice &dipoles=* new matrice(argv[4]);
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
            for(int j=0;j<3;j++) r[j]=dipoles(i,j);
            for(int j=3;j<6;j++) q[j-3]=dipoles(i,j);
            Rs.push_back(r); Qs.push_back(q);
        }

        int newtaille=taille-(geo.getM(geo.nb()-1)).nbTrgs();
        matrice rhs(newtaille, nd);
        assemble_RHS_dipoles( geo, Rs, Qs, rhs,GaussOrder);

        // Saving RHS matrix for dipolar case :
        rhs.SAVE(argv[5]);
    }

	/*********************************************************************************************
    * RK: Computation of RHS for discrete dipolar case: gradient wrt dipoles position and intensity!
    **********************************************************************************************/
    else if(!strcmp(argv[1],"-rhsPOINTgrad")){

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
        int taille=geo.read(argv[2],argv[3]);

        // Loading matrix of dipoles :
        matrice &dipoles=* new matrice(argv[4]);
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
            for(int j=0;j<3;j++) r[j]=dipoles(i,j);
            for(int j=3;j<6;j++) q[j-3]=dipoles(i,j);
            Rs.push_back(r); Qs.push_back(q);
        }

        int newtaille = taille-(geo.getM(geo.nb()-1)).nbTrgs();
        matrice rhs(newtaille, 6*nd); // 6 derivatives! (
        assemble_RHS_dipoles_grad( geo, Rs, Qs, rhs,GaussOrder);

        // Saving RHS matrix for dipolar case :
        rhs.SAVE(argv[5]);
    }

    /*********************************************************************************************
    * Computation of the linear application which maps x (the unknown vector in symmetric system)
    * |----> v (potential at the electrodes)
    **********************************************************************************************/
    else if(!strcmp(argv[1],"-vToEEG"))
    {

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
        int taille = geo.read(argv[2],argv[3]);

        // read the file containing the positions of the EEG patches
        matrice patchesPositions(argv[4]);

        // Assembling matrix from discretization :
        int newtaille = taille-(geo.getM(geo.nb()-1)).nbTrgs();
        matrice vToEEG(patchesPositions.nlin(),newtaille);
        vToEEG.set(0.0);
        assemble_vToEEG( geo, vToEEG, patchesPositions );
        //vToEEG is the linear application which maps x |----> v

        // Saving vToEEG matrix :
        vToEEG.SAVE(argv[5]);
    }

    /*********************************************************************************************
    * Computation of the linear application which maps x (the unknown vector in symmetric system)
    * |----> bFerguson (contrib to MEG response)
    **********************************************************************************************/
    else if(!strcmp(argv[1],"-vToMEG"))
    {

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
        int taille = geo.read(argv[2],argv[3]);

        // Load positions and orientations of sensors  :
        Sensors fileDescription(argv[4]);
        matrice& squidsPositions = fileDescription.getSensorsPositions();
        matrice& squidsOrientations = fileDescription.getSensorsOrientations();

        // Assembling matrix from discretization :
        int newtaille = taille-(geo.getM(geo.nb()-1)).nbTrgs();
        matrice xToMEGrespCont(squidsPositions.nlin(),newtaille);
        xToMEGrespCont.set(0.0);
        assemble_vToMEG( geo, xToMEGrespCont, squidsPositions, squidsOrientations);

        // Saving xToMEGrespCont matrix :
        xToMEGrespCont.SAVE(argv[5]); // if outfile is specified
    }


    /*********************************************************************************************
    * Computation of the linear application which maps x (the unknown vector in symmetric system)
    * |----> binf (contrib to MEG response)
    **********************************************************************************************/
    else if(!strcmp(argv[1],"-sToMEG"))
    {

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
        mesh_sources.load(argv[2],false); // Load mesh without crashing when the surface is not closed

        // Load positions and orientations of sensors  :
        Sensors fileDescription(argv[3]);
        matrice squidsPositions = fileDescription.getSensorsPositions();
        matrice squidsOrientations = fileDescription.getSensorsOrientations();

        // Assembling matrix from discretization :
        int nVertices=mesh_sources.nbPts();
        matrice sToMEGrespCont(squidsPositions.nlin(),nVertices);
        sToMEGrespCont.set(0.0);
        assemble_sToMEG( mesh_sources, sToMEGrespCont, squidsPositions, squidsOrientations);

        // Saving sToMEGrespCont matrix :
        sToMEGrespCont.SAVE(argv[4]);
    }

    /*********************************************************************************************
    * Computation of the discrete linear application which maps x (the unknown vector in a symmetric system)
    * |----> binf (contrib to MEG response)
    **********************************************************************************************/
    // arguments are the positions and orientations of the squids,
    // the position and orientations of the sources and the output name.

    else if(!strcmp(argv[1],"-sToMEG_point"))
    {

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
        size_t nVertices = dipoles.nlin();

        // Load positions and orientations of sensors  :
        Sensors fileDescription(argv[3]);
        matrice squidsPositions = fileDescription.getSensorsPositions();
        matrice squidsOrientations = fileDescription.getSensorsOrientations();

        matrice sToMEGrespCont(squidsPositions.nlin(),nVertices);
        sToMEGrespCont.set(0.0);

        assemble_sToMEG_point( dipoles, sToMEGrespCont, squidsPositions, squidsOrientations);

        sToMEGrespCont.SAVE(argv[4]);
    }
    else cerr << "unknown argument: " << argv[1] << endl;

    // Stop Chrono
    C.stop();
    C.dispEllapsed();
}

void getOutputFilepath(char* ref_filepath, char* output_filename, char* path)
{
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

void getHelp(char** argv)
{
    cout << argv[0] <<" [-option] [filepaths...]" << endl << endl;

    cout << "-option :" << endl;
    cout << "   -LHS :   Compute the LHS from BEM symmetric formulation. " << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            geometry file (.geom), conductivity file (.cond),"  << endl;
    cout << "            name of the output file of LHS matrix (.bin or .txt)" << endl << endl;

    cout << "   -RHS :   Compute the RHS from BEM symmetric formulation. " << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            geometry file (.geom), conductivity file (.cond), " << endl;
    cout << "            mesh file for distributed sources" << endl;
    cout << "            (.tri or .vtk. or .geo), name of the output file of RHS" << endl;
    cout << "            matrix (.bin or .txt)" << endl << endl;

    cout << "   -RHS2 :  Compute in another way (with precomputation of certain" << endl;
    cout << "            blocks in order to speed-up the computation of the " << endl;
    cout << "            diagonal blocks) the RHS from BEM symmetric " << endl;
    cout << "            formulation. Filepaths are in order :" << endl;
    cout << "            geometry file (.geom), conductivity file (.cond), " << endl;
    cout << "            mesh file for distributed sources" << endl;
    cout << "            (.tri or .vtk. or .geo), name of the output file of RHS" << endl;
    cout << "            matrix (.bin or .txt)" << endl << endl;

    cout << "   -rhsPOINT :   Compute the RHS for discrete dipolar case. " << endl;
    cout << "            Filepaths are in order :" << endl;
    cout << "            geometry file (.geom), conductivity file (.cond), " << endl;
    cout << "            file which contains the matrix of " << endl;
    cout << "            dipoles, name of the output file of RHS matrix (.bin or " << endl;
    cout << "            .txt)" << endl << endl;

    cout << "   -vToEEG :   Compute the linear application which maps x " << endl;
    cout << "            (the unknown vector in symmetric system) |----> v (potential"<< endl;
    cout << "            at the electrodes). Filepaths are in order :" << endl;
    cout << "            geometry file (.geom), conductivity file (.cond), " << endl;
    cout << "            file containing the positions of " << endl;
    cout << "            the EEG patches (.patches), name of the output file of " << endl;
    cout << "            xToEEG matrix (.bin or .txt)" << endl << endl;

    cout << "   -vToMEG :   Compute the linear application which maps" << endl;
    cout << "            x (the unknown vector in symmetric system) |----> bFerguson"  << endl;
    cout << "            (contrib to MEG response). Filepaths are in order :" << endl;
    cout << "            geometry file (.geom), conductivity file (.cond), " << endl;
    cout << "            file containing the positions and" << endl;
    cout << "            the orientations of the MEG captors (.squids), name of " << endl;
    cout << "            the output file of xToMEG matrix (.bin or .txt)" << endl << endl;

    cout << "   -sToMEG :   Compute the linear application which maps" << endl;
    cout << "            x (the unknown vector in symmetric system) |----> binf " << endl;
    cout << "            (contrib to MEG response). Filepaths are in order :" << endl;
    cout << "            mesh file for distributed sources (.tri or .vtk. or " << endl;
    cout << "            .geo), file containing the positions and the orientations" << endl;
    cout << "            of the MEG sensors (.squids), name of the output file of" << endl;
    cout << "            sToMEG matrix (.bin or .txt)" << endl << endl;

    cout << "   -sToMEG_point :   Compute the linear application which maps"  << endl;
    cout << "            x (the unknown vector in symmetric system) |----> binf " << endl;
    cout << "            (contrib to MEG response) in dipolar case. Filepaths are in order :" << endl;
    cout << "            mesh file for distributed sources (.tri or .vtk. or " << endl;
    cout << "            .geo), file containing the matrix of dipoles, name of the output " << endl;
    cout << "            file of sToMEG matrix (.bin or .txt)" << endl << endl;

    exit(0);
}

