/*
OpenMEEG

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

#include "assemble.h"
#include "gain.h"
#include "mesh3.h"

using namespace OpenMEEG;

void getHelp(char** argv);

void PotAtInfinity(const Geometry& , const Matrix&, const Matrix&, Matrix&);

int main(int argc, char** argv)
{
    print_version(argv[0]);

    if(argc < 2)
    {
        std::cerr << "Not enough arguments \nPlease try \"" << argv[0] << " -h\" or \"" << argv[0] << " --help \" \n" << std::endl;
        return 0;
    }
    if ((!strcmp(argv[1],"-h")) | (!strcmp(argv[1],"--help"))) {
        getHelp(argv);
        return 0;
    }
    Geometry geo;
    geo.read(argv[1],argv[2]);

    Matrix dipoles;
    dipoles.load(argv[3]);

    Matrix SourceMatrix;
    SourceMatrix.load(argv[4]);

    SymMatrix HeadMatInv;
    HeadMatInv.load(argv[5]);

    Matrix Points;
    std::string extension = getNameExtension(argv[6]); // We check whether it is a mesh file or not
    std::transform(extension.begin(), extension.end(), extension.begin(), (int(*)(int))std::tolower);
    if ((extension==std::string("vtk")) || (extension==std::string("tri")) || (extension==std::string("bnd")) || (extension==std::string("mesh")) || (extension==std::string("off"))){
        Mesh innermesh;
        innermesh.load(argv[6]);
        Matrix pts(innermesh.nbPts(),3);
        for (int i=0;i<innermesh.nbPts();i++){
            for (int j=0;j<3;j++){
                pts(i,j)=(innermesh.getPt(i))(j);
            }
        }
        Points=pts;
    }
    else{ //else we suppose it is a Txt file, with coordinates
        std::cout << "points" << std::endl;
        Points.load(argv[6]);
        if (Points.ncol()!=3){
            std::cerr << "Not a correct file with points coordinates " << argv[6] << std::endl;
            exit(1);
        }
    }

    Surf2VolMat surf2volmat(geo,Points);

    Matrix EEGGainMatrix;
    EEGGainMatrix=surf2volmat*HeadMatInv(0,surf2volmat.ncol()-1,0,HeadMatInv.ncol()-1);
    EEGGainMatrix=EEGGainMatrix*SourceMatrix;

    Matrix Vinfinite(EEGGainMatrix.nlin(),EEGGainMatrix.ncol()); // We must substract the contribution of the dipole at infinity
    PotAtInfinity(geo,Points,dipoles,Vinfinite);
    EEGGainMatrix+=Vinfinite;
    EEGGainMatrix.save(argv[7]);
}

void getHelp(char** argv)
{
    std::cout << "Testing the SurfToVol : \nCompute the potential at points located in the first volume, with coordinates described in the mesh file." << std::endl;

    std::cout << argv[0] << " [filepaths...]" << std::endl;
    std::cout << "Arguments :"                << std::endl;
    std::cout << "               geometry file (.geom)"     << std::endl;
    std::cout << "               conductivity file (.cond)" << std::endl;
    std::cout << "               dipoles file (.dip)"       << std::endl;
    std::cout << "               SourceMat"                 << std::endl;
    std::cout << "               HeadMatInv"                << std::endl ;
    std::cout << "               innermesh or file with the points coordinates on which you want the potential" << std::endl;
    std::cout << "               output Potential" << std::endl;

    exit(0);
}

void PotAtInfinity(const Geometry& geo, const Matrix& pts, const Matrix& dipoles, Matrix &Vinfinite) // Vinf(r)=1/(4*pi*sigma)*(r-r0).q/(||r-r0||^3)
{
    double K = 1/(4*M_PI);

    for (int iDIP=0;iDIP<dipoles.nlin();iDIP++){
        Vect3 r0;
        r0(0)=dipoles(iDIP,0);
        r0(1)=dipoles(iDIP,1);
        r0(2)=dipoles(iDIP,2);
        Vect3 q;
        q(0)=dipoles(iDIP,3);
        q(1)=dipoles(iDIP,4);
        q(2)=dipoles(iDIP,5);
        for (int iPTS=0;iPTS<pts.nlin();iPTS++){
            Vect3 r;
            r(0)=pts(iPTS,0);
            r(1)=pts(iPTS,1);
            r(2)=pts(iPTS,2);
            Vinfinite(iPTS,iDIP)=K*1.0/geo.sigma_in(0)*((r-r0)*q)/(pow((r-r0).norm(),3));
        }
    }
}
