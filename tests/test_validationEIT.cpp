/*
OpenMEEG

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

// This program is a test used to validate the EIT.
// see 'Boundary Element Method for Electrical Impedance Tomography' for more details.
// We check if q.Grad(Vj(r0)) equals Vf(ri)-Vf(re)
// Vf solution of classical forward problem for a dipole of momentum q at r0 whereas
// Vj is the solution of the EIT problem for injections at ri and re

#include "assemble.h"
#include "danielsson.h"

#include "commandline.h"

using namespace OpenMEEG;

void
getHelp(const char* argv[]) {
    std::cout << "Testing the EIT : \n using the Helmholtz reciprocity principle." << std::endl
              << argv[0] << " [filepaths...]" << std::endl
              << "Arguments :"                << std::endl
              << "               geometry file (.geom)"     << std::endl
              << "               conductivity file (.cond)" << std::endl
              << "               dipoles file (.dip)"       << std::endl
              << "               SourceMat"                 << std::endl
              << "               HeadMatInv"                << std::endl
              << "               output q.grad(Vj)         " << std::endl
              << "               output Vf(ri) - Vf(re)      " << std::endl;
    exit(0);
}

Vector
VR(const Geometry& geo,const Matrix& points,const SymMatrix& HeadMatInv,const Matrix& rhsEIT) {
    const Matrix& matdx         = Surf2VolMat(geo,points);
    const Matrix& EEGGainMatrix = matdx*HeadMatInv;
    return (EEGGainMatrix*rhsEIT).getcol(0);
}

int
main(const int argc,const char* argv[]) {

    print_version(argv[0]);

    if (argc<2) {
        std::cerr << "Not enough arguments \nPlease try \"" << argv[0];
        std::cerr << " -h\" or \"" << argv[0] << " --help \" \n" << std::endl;
        return 0;
    }

    if ((!strcmp(argv[1],"-h")) || (!strcmp(argv[1],"--help"))) {
        getHelp(argv);
        return 0;
    }

    Geometry geo(argv[1],argv[2]);

    const unsigned totalsize  = geo.nb_parameters();
    const unsigned sourcesize = geo.outermost_interface().nb_triangles();
    const unsigned newsize    = totalsize-sourcesize;

    const Matrix    dipoles(argv[3]);
    const Matrix    SourceMatrix(argv[4]);
    const SymMatrix HeadMatInv(argv[5]);

    const unsigned ndip = dipoles.nlin();

    // We choose two electrodes on which we inject the currents

    constexpr unsigned nelec = 2;
    const double elec_pos[nelec][3] = {
        { 0.886556, 0.278249,  -0.166667 },
        { 0.547922, -0.269672, -0.719889 }
    };
    Matrix electrodes_positions(2,3);
    for (Index i=0; i<nelec; ++i)
        for (Index j=0; j<3; ++j)
            electrodes_positions(i,j) = elec_pos[i][j];

    Sensors electrodes(electrodes_positions,geo); // set nelec electrode positions

    // We want the potential on the external surface

    HeadMatInv.info();
    SourceMatrix.info();
    const Matrix& PotExt = HeadMatInv*SourceMatrix;
    PotExt.info();

    // Find the triangle closest to the electrodes

    Vect3 current_position; // buffer for electrodes positions
    Vect3 current_alphas;
    Triangle current_nearest_triangle; // buffer for closest triangle to electrode
    SparseMatrix matH2E(nelec,newsize); // Matrices Head2Electrodes

    for (unsigned ielec=0; ielec<nelec; ++ielec) {
        for (unsigned k=0; k<3; ++k)
            current_position(k) = electrodes_positions(ielec,k);
        const auto& res = dist_point_interface(current_position,geo.outermost_interface(),current_alphas);
        const Triangle& current_nearest_triangle = std::get<1>(res); // Closest triangle to electrode
        for (unsigned k=0; k<3; ++k)
            matH2E(ielec,current_nearest_triangle.vertex(k).index()) = current_alphas(k);
    }

    // Potential at the electrodes positions

    const Vector& VRi = (matH2E*PotExt).getlin(0);
    const Vector& VRe = (matH2E*PotExt).getlin(1);

    Matrix diffVf(1,ndip);
    diffVf.setlin(0,VRi-VRe);
    diffVf.save(argv[7]);

    constexpr double dirac = 1.0; // The injection current
    Matrix injection(nelec,1);
    injection(0,0) =  dirac;
    injection(1,0) = -dirac;

    constexpr unsigned gauss_order = 3;
    const Matrix& EITsource = EITSourceMat(geo,electrodes,gauss_order);
    const Matrix& rhsEIT    = EITsource*injection;

    // Surf2Vol
    Matrix points = dipoles.submat(0,ndip,0,3); // extract positions
    const Vector& VR0 = VR(geo,points,HeadMatInv,rhsEIT);

    constexpr double delta = 1e-6;
    Matrix gradVj(ndip,3);
    for (unsigned i=0; i<3; ++i) {
        points.setcol(i,points.getcol(i)+delta);
        const Vector& VRd = VR(geo,points,HeadMatInv,rhsEIT);
        gradVj.setcol(i,(VRd-VR0)/delta);
        points.setcol(i,dipoles.getcol(i));
    }
	
    Matrix qgradVj(1,ndip);
    for (unsigned i=0; i<ndip; ++i)
        qgradVj(0,i) = dipoles(i,3)*gradVj(i,0)+dipoles(i,4)*gradVj(i,1)+dipoles(i,5)*gradVj(i,2);
    qgradVj.save(argv[6]);

    return 0;
}
