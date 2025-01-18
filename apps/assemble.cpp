// Project Name: OpenMEEG (http://openmeeg.github.io)
// © INRIA and ENPC under the French open source license CeCILL-B.
// See full copyright notice in the file LICENSE.txt
// If you make a copy of this file, you must either:
// - provide also LICENSE.txt and modify this header to refer to it.
// - replace this header by the LICENSE.txt content.

#include <fstream>
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

void help(const char* cmd_name);

static inline bool
check_option_is_in_list(const char* option,const std::vector<std::string>& optlist) {
    for (const std::string& opt : optlist)
        if (opt==option)
            return true;
    return false;
}

int main(int argc, char** argv)
{
    print_version(argv[0]);

    const CommandLine cmd(argc,argv,"Compute various head matrices [options] geometry");
    const bool use_old_ordering = cmd.option("-old-ordering", false,"Using old ordering i.e using (V1, p1, V2, p2, V3) instead of (V1, V2, V3, p1, p2)");

    if (argc<2 || cmd.help_mode()) {
        help(argv[0]);
        return 0;
    }

    cmd.print();

    constexpr char geomfileopt[]       = "geometry file";
    constexpr char condfileopt[]       = "conductivity file";
    constexpr char outputfileopt[]     = "output file";
    constexpr char monopolefileopt[]   = "monopoles file";
    constexpr char dipolefileopt[]     = "dipoles file";
    constexpr char electrodesfileopt[] = "electrodes positions file";
    constexpr char squidsfileopt[]     = "squids file";
    constexpr char sourcemeshfileopt[] = "'mesh sources' file";
    constexpr char pointsfileopt[]     = "point positions file";

    // Start Chrono

    const auto start_time = std::chrono::system_clock::now();

    // Computation of Head Matrix for BEM Symmetric formulation

    unsigned num_options = 0;


    const auto& HMparms = { geomfileopt, condfileopt, outputfileopt };
    if (char** opt_parms = cmd.option({ "-HeadMat", "-HM", "-hm" },HMparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        const Geometry geo(opt_parms[1],opt_parms[2],use_old_ordering);

        if (!geo.selfCheck()) // Check for intersecting meshes
            exit(1);

        const SymMatrix& HM = HeadMat(geo);
        HM.save(opt_parms[3]);
    }

    const auto& CMparms = { geomfileopt, condfileopt, "sensors file", "domain name", outputfileopt };
    if (char** opt_parms = cmd.option({ "-CorticalMat", "-CM", "-cm" },CMparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        //  Computation of Cortical Matrix for BEM Symmetric formulation

        double alpha = -1.0;
        double beta  = -1.0;
        double gamma = -1.0;

        std::string filename;

        switch (cmd.num_args(opt_parms)) {
            case 6: { // case gamma or filename
                std::stringstream ss(opt_parms[6]);
                if (!(ss >> gamma))
                    filename = opt_parms[6];
                break;
            }
            case 7: { // case alpha+beta or gamma+filename
                std::stringstream ss(opt_parms[6]);
                if (!(ss >> alpha))
                    throw std::runtime_error("given parameter is not a number");
                ss.str(opt_parms[7]);
                ss.clear();
                if (!(ss >> beta)) {
                    filename = opt_parms[7];
                    gamma = alpha;
                }
                break;
            }
            case 8: { // case alpha+beta + filename
                std::stringstream ss(opt_parms[6]);
                if (!(ss >> alpha))
                    throw std::runtime_error("given parameter is not a number");
                ss.str(opt_parms[7]);
                ss.clear();
                if (!(ss >> beta))
                    throw std::runtime_error("given parameter is not a number");
                filename = opt_parms[8];
                break;
            }
        }

        const Geometry geo(opt_parms[1],opt_parms[2],use_old_ordering);

        if (!geo.selfCheck()) // Check for intersecting meshes
            exit(1);

        const Sensors electrodes(opt_parms[3]);
        const SparseMatrix& M = Head2EEGMat(geo,electrodes);

        const Matrix& CM = (gamma>0.0) ? CorticalMat2(geo,M,opt_parms[4],gamma,filename) :
                                         CorticalMat(geo,M,opt_parms[4],alpha,beta,filename);
        CM.save(opt_parms[5]);
    }

    const auto& SSMparms = { geomfileopt, condfileopt, sourcemeshfileopt, outputfileopt };
    if (char** opt_parms = cmd.option({"-SurfSourceMat", "-SSM", "-ssm"},SSMparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        // Computation of distributed Surface Source Matrix for BEM Symmetric formulation

        const Geometry geo(opt_parms[1],opt_parms[2],use_old_ordering);
        Mesh mesh_sources(opt_parms[3]);

        const Matrix& ssm = SurfSourceMat(geo,mesh_sources);
        ssm.save(opt_parms[4]);
    }

    const auto& MSMparms = { geomfileopt, condfileopt, monopolefileopt, outputfileopt };
    if (char** opt_parms = cmd.option({"-MonopoleSourceMat", "-MSM", "-msm", "-MonopoleSourceMatNoAdapt", "-MSMNA", "-msmna"},MSMparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        // Computation of RHS for discrete monopolar case

        std::string domain_name = "";
        if (cmd.num_args(opt_parms)==5) {
            domain_name = opt_parms[6];
            std::cout << "Monopoles are considered to be in \"" << domain_name << "\" domain." << std::endl;
        }

        const Geometry geo(opt_parms[1],opt_parms[2],use_old_ordering);
        const Matrix monopoles(opt_parms[3]);
        if (monopoles.ncol()!=4) {
            std::cerr << "Monopoles File Format Error" << std::endl;
            exit(1);
        }

        // Choosing between adaptive integration or not for the RHS

        const char* optname = opt_parms[0];
        const unsigned integration_levels = check_option_is_in_list(optname,{"-MonopoleSourceMatNoAdapt", "-MSMNA", "-msmna"}) ? 0 : 10;

        const Matrix& msm = MonopoleSourceMat(geo,monopoles,Integrator(3,integration_levels,0.001),domain_name);
        msm.save(opt_parms[4]);
    }

    const auto& DSMparms = { geomfileopt, condfileopt, dipolefileopt, outputfileopt };
    if (char** opt_parms = cmd.option({"-DipSourceMat", "-DSM", "-dsm", "-DipSourceMatNoAdapt", "-DSMNA", "-dsmna"},DSMparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        // Computation of RHS for discrete dipolar case

        std::string domain_name = "";
        if (cmd.num_args(opt_parms)==5) {
            domain_name = opt_parms[6];
            std::cout << "Dipoles are considered to be in \"" << domain_name << "\" domain." << std::endl;
        }

        const Geometry geo(opt_parms[1],opt_parms[2],use_old_ordering);
        const Matrix dipoles(opt_parms[3]);
        if (dipoles.ncol()!=6) {
            std::cerr << "Dipoles File Format Error" << std::endl;
            exit(1);
        }

        // Choosing between adaptive integration or not for the RHS

        const char* optname = opt_parms[0];
        const unsigned integration_levels = check_option_is_in_list(optname,{"-DipSourceMatNoAdapt", "-DSMNA", "-dsmna"}) ? 0 : 10;

        const Matrix& dsm = DipSourceMat(geo,dipoles,Integrator(3,integration_levels,0.001),domain_name);
        dsm.save(opt_parms[4]);
    }

    const auto& EITSMparms = { geomfileopt, condfileopt, electrodesfileopt, outputfileopt };
    if (char** opt_parms = cmd.option({"-EITSourceMat", "-EITSM", "-EITsm"},EITSMparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        // Computation of the RHS for EIT

        const Geometry geo(opt_parms[1],opt_parms[2],use_old_ordering);

        const Sensors electrodes(opt_parms[3], geo); // special parameter for EIT electrodes: the interface
        electrodes.info(); // <- just to test that function on the code coverage TODO this is not the place.
        const Matrix& EITsource = EITSourceMat(geo,electrodes);
        EITsource.save(opt_parms[4]);
    }

    const auto& H2EMparms = { geomfileopt, condfileopt, electrodesfileopt, outputfileopt };
    if (char** opt_parms = cmd.option({"-Head2EEGMat", "-H2EM", "-h2em"},H2EMparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        // Computation of the linear application which maps the unknown vector in symmetric system,
        // (i.e. the potential and the normal current on all interfaces)
        // |----> v (potential at the electrodes)

        const Geometry geo(opt_parms[1],opt_parms[2],use_old_ordering);
        const Sensors electrodes(opt_parms[3]);

        // Head2EEG is the linear application which maps x |----> v

        const SparseMatrix& mat = Head2EEGMat(geo,electrodes);
        mat.save(opt_parms[4]);
    }

    const auto& H2ECOGMparms = {
        geomfileopt, condfileopt, "ECoG electrodes positions file", "[name of the interface for EcoG]", outputfileopt
    };
    if (char** opt_parms = cmd.option({"-Head2ECoGMat", "-H2ECogM", "-H2ECOGM", "-h2ecogm"},H2ECOGMparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        // Computation of the linear application which maps the unknown vector in symmetric system,
        // (i.e. the potential and the normal current on all interfaces)
        // |----> v (potential at the ECoG electrodes)

        const Geometry geo(opt_parms[1],opt_parms[2],use_old_ordering);
        const Sensors electrodes(opt_parms[3]);

        // Find the mesh of the ECoG electrodes

        const bool old_cmd_line = (cmd.num_args(opt_parms)==4) ? true : false;
        if (old_cmd_line)
            std::cerr << "Warning: we assume that ECoG electrodes are placed on the inner interface." << std::endl
                      << "This is only valid for nested files. Consider specifying an interface as a name" << std::endl
                      << " right after the electrode position file." << std::endl;

        const Interface& ECoG_layer = (old_cmd_line) ? geo.innermost_interface() : geo.interface(opt_parms[4]);

        // Head2ECoG is the linear application which maps x |----> v

        const SparseMatrix& mat = Head2ECoGMat(geo,electrodes,ECoG_layer);
        mat.save(opt_parms[(old_cmd_line) ? 4 : 5]);
    }

    const auto& H2MMparms = { geomfileopt, condfileopt, squidsfileopt, outputfileopt };
    if (char** opt_parms = cmd.option({"-Head2MEGMat", "-H2MM", "-h2mm"},H2MMparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        // Computation of the linear application which maps the unknown vector in symmetric system,
        // (i.e. the potential and the normal current on all interfaces)
        // |----> bFerguson (contrib to MEG response)

        const Geometry geo(opt_parms[1],opt_parms[2],use_old_ordering);
        const Sensors sensors(opt_parms[3]);

        const Matrix& mat = Head2MEGMat(geo,sensors);
        mat.save(opt_parms[4]); // if outfile is specified
    }

    const auto& SS2MMparms = { sourcemeshfileopt, squidsfileopt, outputfileopt };
    if (char** opt_parms = cmd.option({"-SurfSource2MEGMat", "-SS2MM", "-ss2mm"},SS2MMparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        // Computation of the linear application which maps the distributed source
        // |----> binf (contrib to MEG response)

        const Mesh mesh_sources(opt_parms[1]);
        const Sensors sensors(opt_parms[2]);

        const Matrix& mat = SurfSource2MEGMat(mesh_sources,sensors);
        mat.save(opt_parms[3]);
    }

    const auto& DS2MMparms = { dipolefileopt, squidsfileopt, outputfileopt };
    if (char** opt_parms = cmd.option({"-DipSource2MEGMat", "-DS2MM", "-ds2mm"},DS2MMparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        // Computation of the discrete linear application which maps s (the dipolar source)
        // |----> binf (contrib to MEG response)
        // arguments are the positions and orientations of the squids,
        // the position and orientations of the sources and the output name.

        const Matrix dipoles(opt_parms[1]);
        const Sensors sensors(opt_parms[2]);

        const Matrix& mat = DipSource2MEGMat(dipoles,sensors);
        mat.save(opt_parms[3]);
    }

    const auto& H2IPMparms = { geomfileopt, condfileopt, pointsfileopt, outputfileopt };
    if (char** opt_parms = cmd.option({"-Head2InternalPotMat", "-H2IPM", "-h2ipm"},H2IPMparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        // Computation of the discrete linear application which maps x (the unknown vector in a symmetric system)
        // |----> v, potential at a set of prescribed points within the 3D volume

        const Geometry geo(opt_parms[1],opt_parms[2],use_old_ordering);
        const Matrix points(opt_parms[3]);
        const Matrix& mat = Surf2VolMat(geo,points);
        mat.save(opt_parms[4]);
    }

    const auto& DS2IPMparms = { geomfileopt, condfileopt, dipolefileopt, pointsfileopt, outputfileopt };
    if (char** opt_parms = cmd.option({"-DipSource2InternalPotMat", "-DS2IPM", "-ds2ipm"},DS2IPMparms)) {

        assert_non_conflicting_options(argv[0],++num_options);

        // Computation of the discrete linear application which maps the dipoles
        // |----> Vinf, potential at a set of prescribed points within the volume, in an infinite medium
        //    Vinf(r)=1/(4*pi*sigma)*(r-r0).q/(||r-r0||^3)

        const std::string& domain_name = (cmd.num_args(opt_parms)==6) ? opt_parms[6] : "";
        if (domain_name!="")
            std::cout << "Dipoles are considered to be in \"" << domain_name << "\" domain." << std::endl;

        const Geometry geo(opt_parms[1],opt_parms[2],use_old_ordering);
        const Matrix dipoles(opt_parms[3]);
        const Matrix points(opt_parms[4]);
        const Matrix& mat = DipSource2InternalPotMat(geo,dipoles,points,domain_name);
        mat.save(opt_parms[5]);
    }

    if (num_options==0) {
        std::cerr << "unknown argument: " << argv[1] << std::endl;
        exit(1);
    }

    // Stop Chrono
    const auto end_time = std::chrono::system_clock::now();
    dispEllapsed(end_time-start_time);

    return 0;
}

void help(const char* cmd_name) {
    std::cout << cmd_name <<" [-option] [filepaths...]" << std::endl << std::endl;

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

    std::cout << "   -MonopoleSourceMat, -MSM, -msm:    " << std::endl
              << "      Compute Monopolar Source Matrix for Symmetric BEM (right-hand side of linear system). " << std::endl
              << "            Arguments:" << std::endl
              << "               geometry file (.geom)" << std::endl
              << "               conductivity file (.cond)" << std::endl
              << "               monopole positions and charges" << std::endl
              << "               output matrix" << std::endl
              << "               (Optional) domain name where lie all monopoles." << std::endl << std::endl;

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
}
