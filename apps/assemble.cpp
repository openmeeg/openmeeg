/* OpenMEEG

© INRIA and ENPC (contributors: Geoffray ADDE, Maureen CLERC, Alexandre
GRAMFORT, Renaud KERIVEN, Jan KYBIC, Perrine LANDREAU, Théodore PAPADOPOULO,
Emmanuel OLIVI
Maureen.Clerc.AT.inria.fr, keriven.AT.certis.enpc.fr,
kybic.AT.fel.cvut.cz, papadop.AT.inria.fr)

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

#include <fstream>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <mesh.h>
#include <integrator.h>
#include <om_utils.h>
#include <commandline.h>
#include <assemble.h>
#include <sensors.h>
#include <geometry.h>

using namespace OpenMEEG;

void getHelp(char** argv);

bool option(const int argc,char** argv,const Strings& options,const Strings& files);

int main(int argc, char** argv)
{
    print_version(argv[0]);

    bool use_old_ordering = false;
    if (argc<2) {
        getHelp(argv);
        return 0;
    } else {
        use_old_ordering = (strcmp(argv[argc-1], "-old-ordering") == 0);
        if (use_old_ordering) {
            argc--;
            std::cout << "Using old ordering i.e using (V1, p1, V2, p2, V3) instead of (V1, V2, V3, p1, p2)" << std::endl;
        }
    }

    if (option(argc,argv,{"-h","--help"}, {})) getHelp(argv);

    print_commandline(argc, argv);

    // Start Chrono

    auto start_time = std::chrono::system_clock::now();

    /*********************************************************************************************
    * Computation of Head Matrix for BEM Symmetric formulation
    **********************************************************************************************/

    if (option(argc,argv,{"-HeadMat","-HM", "-hm"},
                {"geometry file", "conductivity file", "output file"}) ) {
        const Geometry geo(argv[2],argv[3],use_old_ordering);

        if (!geo.selfCheck()) // Check for intersecting meshes
            exit(1);

        const SymMatrix& HM = HeadMat(geo);
        HM.save(argv[4]);

    } else if (option(argc,argv,{ "-CorticalMat","-CM","-cm" },
                                { "geometry file","conductivity file","sensors file","domain name","output file" })) {

        //  Computation of Cortical Matrix for BEM Symmetric formulation

        double alpha = -1.0;
        double beta  = -1.0;
        double gamma = -1.0;

        std::string filename;

        switch (argc) {
            case 8: { // case gamma or filename
                std::stringstream ss(argv[7]);
                if (!(ss >> gamma))
                    filename = argv[7];
                break;
            }
            case 9: { // case alpha+beta or gamma+filename
                std::stringstream ss(argv[7]);
                if (!(ss >> alpha))
                    throw std::runtime_error("given parameter is not a number");
                ss.str(argv[8]);
                ss.clear();
                if (!(ss >> beta)) {
                    filename = argv[8];
                    gamma = alpha;
                }
                break;
            }
            case 10: { // case alpha+beta + filename
                std::stringstream ss(argv[7]);
                if (!(ss >> alpha))
                    throw std::runtime_error("given parameter is not a number");
                ss.str(argv[8]);
                ss.clear();
                if (!(ss >> beta))
                    throw std::runtime_error("given parameter is not a number");
                filename = argv[9];
                break;
            }
        }

        const Geometry geo(argv[2],argv[3],use_old_ordering);

        if (!geo.selfCheck()) // Check for intersecting meshes
            exit(1);

        const Sensors electrodes(argv[4]);
        const SparseMatrix& M = Head2EEGMat(geo,electrodes);

        const Matrix& CM = (gamma>0.0) ? CorticalMat2(geo,M,argv[5],gamma,filename) :
                                         CorticalMat(geo,M,argv[5],alpha,beta,filename);
        CM.save(argv[6]);

    } else if (option(argc,argv,{"-SurfSourceMat", "-SSM", "-ssm"},
                     {"geometry file", "conductivity file", "'mesh of sources' file", "output file"})) {

        // Computation of distributed Surface Source Matrix for BEM Symmetric formulation
        
        const Geometry geo(argv[2],argv[3],use_old_ordering);
        Mesh mesh_sources(argv[4]);

        const Matrix& ssm = SurfSourceMat(geo,mesh_sources);
        ssm.save(argv[5]);

    } else if (option(argc,argv,{"-DipSourceMat", "-DSM", "-dsm", "-DipSourceMatNoAdapt", "-DSMNA", "-dsmna"},
                     {"geometry file", "conductivity file", "dipoles file", "output file"}) ) {

        // Computation of RHS for discrete dipolar case
    
        std::string domain_name = "";
        if (argc==7) {
            domain_name = argv[6];
            std::cout << "Dipoles are considered to be in \"" << domain_name << "\" domain." << std::endl;
        }

        const Geometry geo(argv[2],argv[3],use_old_ordering);
        const Matrix dipoles(argv[4]);
        if (dipoles.ncol()!=6) {
            std::cerr << "Dipoles File Format Error" << std::endl;
            exit(1);
        }

        // Choosing between adaptive integration or not for the RHS

        unsigned integration_levels = 10;
        if (option(argc,argv,{"-DipSourceMatNoAdapt", "-DSMNA", "-dsmna"},
                    {"geometry file", "conductivity file", "dipoles file", "output file"})) {
            integration_levels = 0;
        }

        const Matrix& dsm = DipSourceMat(geo,dipoles,Integrator(3,integration_levels,0.001));
        dsm.save(argv[5]);

    } else if (option(argc,argv,{"-EITSourceMat", "-EITSM", "-EITsm"},
                     {"geometry file", "conductivity file", "electrodes positions file", "output file"}) ) {

        // Computation of the RHS for EIT

        const Geometry geo(argv[2],argv[3],use_old_ordering);

        const Sensors electrodes(argv[4], geo); // special parameter for EIT electrodes: the interface
        electrodes.info(); // <- just to test that function on the code coverage TODO this is not the place.
        const Matrix& EITsource = EITSourceMat(geo,electrodes);
        EITsource.save(argv[5]);

    } else if (option(argc,argv,{"-Head2EEGMat", "-H2EM", "-h2em"},
                     {"geometry file", "conductivity file", "electrodes positions file", "output file"}) ) {

        // Computation of the linear application which maps the unknown vector in symmetric system,
        // (i.e. the potential and the normal current on all interfaces)
        // |----> v (potential at the electrodes)

        const Geometry geo(argv[2],argv[3],use_old_ordering);
        const Sensors electrodes(argv[4]);

        // Head2EEG is the linear application which maps x |----> v

        const SparseMatrix& mat = Head2EEGMat(geo,electrodes);
        mat.save(argv[5]);

    } else if (option(argc, argv, {"-Head2ECoGMat", "-H2ECogM", "-H2ECOGM", "-h2ecogm"},
                                {"geometry file", "conductivity file", "ECoG electrodes positions file",
                                 "[name of the interface for EcoG]", "output file"}) ) {

        // Computation of the linear application which maps the unknown vector in symmetric system,
        // (i.e. the potential and the normal current on all interfaces)
        // |----> v (potential at the ECoG electrodes)

        const Geometry geo(argv[2],argv[3],use_old_ordering);
        const Sensors electrodes(argv[4]);

        // Find the mesh of the ECoG electrodes

        bool old_cmd_line = false;
        if (argc==6) {
            std::cerr << "Warning: we assume that ECoG electrodes are placed on the inner interface." << std::endl
                      << "This is only valid for nested files. Consider specifying an interface as a name" << std::endl
                      << " right after the electrode position file." << std::endl;
            old_cmd_line = true;
        }

        const Interface& ECoG_layer = (old_cmd_line) ? geo.innermost_interface() : geo.interface(argv[5]);

        // Head2ECoG is the linear application which maps x |----> v

        const SparseMatrix& mat = Head2ECoGMat(geo,electrodes,ECoG_layer);
        mat.save(argv[(old_cmd_line) ? 5 : 6]);

    } else if (option(argc,argv,{"-Head2MEGMat", "-H2MM", "-h2mm"},
                     {"geometry file", "conductivity file", "squids file", "output file"}) ) {

        // Computation of the linear application which maps the unknown vector in symmetric system,
        // (i.e. the potential and the normal current on all interfaces)
        // |----> bFerguson (contrib to MEG response)

        const Geometry geo(argv[2],argv[3],use_old_ordering);
        const Sensors sensors(argv[4]);

        const Matrix& mat = Head2MEGMat(geo,sensors);
        mat.save(argv[5]); // if outfile is specified

    } else if (option(argc,argv,{"-SurfSource2MEGMat", "-SS2MM", "-ss2mm"},
                     {"'mesh sources' file", "squids file", "output file"}) ) {

        // Computation of the linear application which maps the distributed source
        // |----> binf (contrib to MEG response)

        const Mesh mesh_sources(argv[2]);
        const Sensors sensors(argv[3]);

        const Matrix& mat = SurfSource2MEGMat(mesh_sources,sensors);
        mat.save(argv[4]);

    } else if (option(argc,argv,{"-DipSource2MEGMat", "-DS2MM", "-ds2mm"},
                     {"dipoles file", "squids file", "output file"}) ) {

        // Computation of the discrete linear application which maps s (the dipolar source)
        // |----> binf (contrib to MEG response)
        // arguments are the positions and orientations of the squids,
        // the position and orientations of the sources and the output name.

        const Matrix dipoles(argv[2]);
        const Sensors sensors(argv[3]);

        const Matrix& mat = DipSource2MEGMat(dipoles,sensors);
        mat.save(argv[4]);

    } else if (option(argc,argv,{"-Head2InternalPotMat", "-H2IPM", "-h2ipm"},
                      {"geometry file", "conductivity file", "point positions file", "output file"}) ) {

        // Computation of the discrete linear application which maps x (the unknown vector in a symmetric system)
        // |----> v, potential at a set of prescribed points within the 3D volume

        const Geometry geo(argv[2],argv[3],use_old_ordering);
        const Matrix points(argv[4]);
        const Matrix& mat = Surf2VolMat(geo,points);
        mat.save(argv[5]);
    }

    else if (option(argc,argv,{"-DipSource2InternalPotMat", "-DS2IPM", "-ds2ipm"},
                    {"geometry file", "conductivity file", "dipole file", "point positions file", "output file"})) {

        // Computation of the discrete linear application which maps the dipoles
        // |----> Vinf, potential at a set of prescribed points within the volume, in an infinite medium
        //    Vinf(r)=1/(4*pi*sigma)*(r-r0).q/(||r-r0||^3)

        std::string domain_name = "";
        if (argc==9) {
            domain_name = argv[7];
            std::cout << "Dipoles are considered to be in \"" << domain_name << "\" domain." << std::endl;
        }
        const Geometry geo(argv[2],argv[3],use_old_ordering);
        const Matrix dipoles(argv[4]);
        const Matrix points(argv[5]);
        const Matrix& mat = DipSource2InternalPotMat(geo,dipoles,points,domain_name);
        mat.save(argv[6]);
    }
    else {
        std::cerr << "unknown argument: " << argv[1] << std::endl;
        exit(1);
    }

    // Stop Chrono
    auto end_time = std::chrono::system_clock::now();
    dispEllapsed(end_time-start_time);

    return 0;
}

bool option(const int argc, char ** argv, const Strings& options, const Strings& files) {
    for ( auto s: options) {
        if (argv[1] == s) {
            unsigned minimum_nparms = files.size();
            for (auto f: files)
                if (f[0]=='[')
                    --minimum_nparms;
            if (argc-2 < minimum_nparms) {
                std::cout << "\'om_assemble\' option \'" << argv[1] << "\' expects " << files.size() << " arguments (";
                for ( auto f: files) {
                    std::cout << f;
                    if (f != files[files.size() - 1]){
                        std::cout << ", ";
                    }
                }
                std::cout << ") and you gave only " << argc-2 << " arguments." << std::endl;
                exit(1);
            }
            return true;
        }
    }
    return false;
}

void getHelp(char** argv) {
    std::cout << argv[0] <<" [-option] [filepaths...]" << std::endl << std::endl;

    std::cout << "option:" << std::endl
              << "   -HeadMat, -HM, -hm:   " << std::endl
              << "       Compute Head Matrix for Symmetric BEM (left-hand side of linear system)." << std::endl
              << "             Arguments:" << std::endl
              << "               geometry file (.geom)" << std::endl
              << "               conductivity file (.cond)" << std::endl
              << "               output matrix" << std::endl << std::endl;

    std::cout << "   -CorticalMat, -CM, -cm:   " << std::endl
              << "       Compute Cortical Matrix for Symmetric BEM (left-hand side of linear system)." << std::endl
              << "       Comment on optional parameters:" << std::endl
              << "       Giving two (or zero) numeric optional parameters => CorticalMat will try to use (or estimate) alpha/beta." << std::endl
              << "       Giving one numeric optional parameters => CorticalMat2 will use gamma." << std::endl
              << "       Giving a filename (a string), one can save time saving the intermediate matrix in all cases (useful when trying multiple values)." << std::endl
              << "             Arguments:" << std::endl
              << "               geometry file (.geom)" << std::endl
              << "               conductivity file (.cond)" << std::endl
              << "               file containing the positions of EEG electrodes (.patches)" << std::endl
              << "               domain name (containing the sources)" << std::endl
              << "               output matrix" << std::endl
              << "               [optional parameter alpha or gamma or filename]" << std::endl
              << "               [optional parameter beta or filename]" << std::endl
              << "               [optional filename]" << std::endl << std::endl;

    std::cout << "   -SurfSourceMat, -SSM, -ssm:   " << std::endl
              << "       Compute Surfacic Source Matrix for Symmetric BEM (right-hand side of linear system). " << std::endl
              << "            Arguments:" << std::endl
              << "               geometry file (.geom)" << std::endl
              << "               conductivity file (.cond)" << std::endl
              << "               mesh of sources (.tri .vtk .mesh .bnd)" << std::endl
              << "               output matrix" << std::endl << std::endl;

    std::cout << "   -DipSourceMat, -DSM, -dsm:    " << std::endl
              << "      Compute Dipolar Source Matrix for Symmetric BEM (right-hand side of linear system). " << std::endl
              << "            Arguments:" << std::endl
              << "               geometry file (.geom)" << std::endl
              << "               conductivity file (.cond)" << std::endl
              << "               dipoles positions and orientations" << std::endl
              << "               output matrix" << std::endl
              << "               (Optional) domain name where lie all dipoles." << std::endl << std::endl;

    std::cout << "   -EITSourceMat, -EITSM -EITsm: " << std::endl
              << "       Compute the EIT Source Matrix from an injected current (right-hand side of linear system). " << std::endl
              << "            Arguments:" << std::endl
              << "               geometry file (.geom)" << std::endl
              << "               conductivity file (.cond)" << std::endl
              << "               file containing the positions of EIT electrodes (.patches)" << std::endl
              << "               output EITSourceOp" << std::endl;

    std::cout << "   -Head2EEGMat, -H2EM, -h2em: " << std::endl
              << "        Compute the linear application which maps the potential" << std::endl
              << "        on the scalp to the EEG electrodes"  << std::endl
              << "            Arguments:" << std::endl
              << "               geometry file (.geom)" << std::endl
              << "               conductivity file (.cond)" << std::endl
              << "               file containing the positions of EEG electrodes (.patches)" << std::endl
              << "               output matrix" << std::endl << std::endl;

    std::cout << "   -Head2ECoGMat, -H2ECogM, -h2ecogm, -H2ECOGM: " << std::endl
              << "        Compute the linear application which maps the potential" << std::endl
              << "        on the scalp to the ECoG electrodes"  << std::endl
              << "            Arguments:" << std::endl
              << "               geometry file (.geom)" << std::endl
              << "               conductivity file (.cond)" << std::endl
              << "               file containing the positions of ECoG electrodes (.patches)" << std::endl
              << "               name of the interface on which to project the electrodes (\"name\")" << std::endl
              << "               output matrix" << std::endl << std::endl;

    std::cout << "   -Head2MEGMat, -H2MM, -h2mm: " << std::endl
              << "        Compute the linear application which maps the potential" << std::endl
              << "        on the scalp to the MEG sensors"  << std::endl
              << "            Arguments:" << std::endl
              << "               geometry file (.geom)" << std::endl
              << "               conductivity file (.cond)" << std::endl
              << "               file containing the positions and orientations of the MEG sensors (.squids)" << std::endl
              << "               output matrix" << std::endl << std::endl;

    std::cout << "   -SurfSource2MEGMat, -SS2MM, -ss2mm: " << std::endl
              << "        Compute the linear application which maps the " << std::endl
              << "        distributed source  to the MEG sensors" << std::endl
              << "            Arguments:" << std::endl
              << "               mesh file for distributed sources (.tri .vtk .mesh .bnd)" << std::endl
              << "               positions and orientations of the MEG sensors (.squids)" << std::endl
              << "               output matrix" << std::endl << std::endl;

    std::cout << "   -DipSource2MEGMat, -DS2MM, -ds2mm:  " << std::endl
              << "        Compute the linear application which maps the current dipoles" << std::endl
              << "        to the MEG sensors" << std::endl
              << "            Arguments:" << std::endl
              << "               dipoles positions and orientations" << std::endl
              << "               positions and orientations of the MEG sensors (.squids)" << std::endl
              << "               output matrix" << std::endl << std::endl;

    std::cout << "   -Head2InternalPotMat, -H2IPM -h2ipm:  " << std::endl
              << "        Compute the linear transformation which maps the surface potential" << std::endl
              << "        and normal current to the value of the internal potential at a set of points within a volume" << std::endl
              << "            Arguments:" << std::endl
              << "               geometry file (.geom)" << std::endl
              << "               conductivity file (.cond)" << std::endl
              << "               a mesh file or a file with point positions at which to evaluate the potential" << std::endl
              << "               output matrix" << std::endl << std::endl;

    std::cout << "   -DipSource2InternalPotMat, -DS2IPM -ds2ipm:   " << std::endl
              << "        Compute the linear transformation  which maps the current dipoles" << std::endl
              << "        to the value of the infinite potential at a set of points within a volume" << std::endl
              << "            Arguments:" << std::endl
              << "               geometry file (.geom)" << std::endl
              << "               conductivity file (.cond)" << std::endl
              << "               dipoles positions and orientations" << std::endl
              << "               a mesh file or a file with point positions at which to evaluate the potential" << std::endl
              << "               output matrix" << std::endl
              << "               (Optional) domain name where lie all dipoles." << std::endl << std::endl;

    exit(0);
}
