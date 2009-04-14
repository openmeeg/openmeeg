#include "assemble.h"
#include "gain.h"
#include "mesh3.h"

#define SAVETXT

#ifdef SAVEBIN
    #define SAVE saveBin
#else
    #define SAVE saveTxt
#endif

using namespace OpenMEEG;

void getHelp(char** argv);

void PotAtInfinity(const Geometry& , const Matrix&, const Matrix&, Matrix&);

Matrix GetPoints();

int main(int argc, char** argv)
{
    if ((!strcmp(argv[1],"-h")) | (!strcmp(argv[1],"--help"))) {
        getHelp(argv);
        return 0;
    }
    Geometry geo;
    geo.read(argv[1],argv[2]);

    Matrix dipoles;
    dipoles.loadTxt(argv[3]);

    Matrix SourceMatrix;
    SourceMatrix.loadBin(argv[4]);

    SymMatrix HeadMatInv;
    HeadMatInv.loadBin(argv[5]);

    Matrix Points;
    if (argc>7){// if a mesh file is given, we compute the potential at these points
        Mesh innermesh;
        innermesh.load(argv[7]);
        Matrix pts(innermesh.nbPts(),3);
        for (int i=0;i<innermesh.nbPts();i++){
            for (int j=0;j<3;j++){
                pts(i,j)=(innermesh.getPt(i))(j);
            }
        }
        Points=pts;
    }
    else{ // else we load 15 points
        Points=GetPoints();
    }

    Surf2VolMat surf2volmat(geo,Points);

    Matrix EEGGainMatrix;
    EEGGainMatrix=surf2volmat*HeadMatInv(0,surf2volmat.ncol()-1,0,HeadMatInv.ncol()-1);
    EEGGainMatrix=EEGGainMatrix*SourceMatrix;

    Matrix Vinfinite(EEGGainMatrix.nlin(),EEGGainMatrix.ncol()); // We must substract the contribution of the dipole at infinity
    PotAtInfinity(geo,Points,dipoles,Vinfinite);
    EEGGainMatrix+=Vinfinite;
    EEGGainMatrix.SAVE(argv[6]);
}

void getHelp(char** argv)
{
    std::cout << "Testing the SurfToVol : \nCompute the potential at points located in the first volume, with coordinates described in the mesh file if given or at 10  inner points." << std::endl;

    std::cout << argv[0] << " [filepaths...]" << std::endl;
    std::cout << "Arguments :"                << std::endl;
    std::cout << "               geometry file (.geom)"     << std::endl;
    std::cout << "               conductivity file (.cond)" << std::endl;
    std::cout << "               dipoles file (.dip)"       << std::endl;
    std::cout << "               SourceMat"                 << std::endl;
    std::cout << "               HeadMatInv"                << std::endl ;
    std::cout << "               output Potential" << std::endl;
    std::cout << "   (optional)  innermesh on which you want the potential" << std::endl;
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
        std:: cout << r0(0) << "," << r0(1) << "," << r0(2) << std::endl;
        for (int iPTS=0;iPTS<pts.nlin();iPTS++){
            Vect3 r;
            r(0)=pts(iPTS,0);
            r(1)=pts(iPTS,1);
            r(2)=pts(iPTS,2);
            Vinfinite(iPTS,iDIP)=K*1.0/geo.sigma_in(0)*((r-r0)*q)/(pow((r-r0).norm(),3));
        }
    }
}

Matrix GetPoints() // We set up 15 points in the inner volume at which we have already computed the analytical solution.
{
    Matrix pts(15,3);
    pts(0,0)=0.08005925;	 pts(0,1)=0.00234002;	 pts(0,2)=0.05987521;	 //r=0.100000 theta=0.029220 phi=0.928854
    pts(1,0)=0.06993597;	 pts(1,1)=0.06262861;	 pts(1,2)=0.17659733;	 //r=0.200000 theta=0.730331 phi=0.488609
    pts(2,0)=0.03936227;	 pts(2,1)=0.02570602;	 pts(2,2)=0.19439602;	 //r=0.200000 theta=0.578525 phi=0.237284
    pts(3,0)=0.28705616;	 pts(3,1)=0.14181029;	 pts(3,2)=0.22268499;	 //r=0.390000 theta=0.458849 phi=0.963089
    pts(4,0)=0.17010865;	 pts(4,1)=0.10354833;	 pts(4,2)=0.34690170;	 //r=0.400000 theta=0.546806 phi=0.521136
    pts(5,0)=0.18284564;	 pts(5,1)=0.04311972;	 pts(5,2)=0.35314043;	 //r=0.400000 theta=0.231594 phi=0.488898
    pts(6,0)=0.22937753;	 pts(6,1)=0.16516472;	 pts(6,2)=0.35015220;	 //r=0.450000 theta=0.624060 phi=0.679136
    pts(7,0)=0.18232037;	 pts(7,1)=0.07612181;	 pts(7,2)=0.51328818;	 //r=0.550000 theta=0.395515 phi=0.367437
    pts(8,0)=0.01245938;	 pts(8,1)=0.01890081;	 pts(8,2)=0.59957278;	 //r=0.600000 theta=0.987982 phi=0.037739
    pts(9,0)=0.30069379;	 pts(9,1)=0.36759053;	 pts(9,2)=0.36668849;	 //r=0.600000 theta=0.885168 phi=0.913287
    pts(10,0)=0.05445562;	 pts(10,1)=0.05564316;	 pts(10,2)=0.78615420;	 //r=0.790000 theta=0.796184 phi=0.098712
    pts(11,0)=0.26066653;	 pts(11,1)=0.06986545;	 pts(11,2)=0.77432020;	 //r=0.820000 theta=0.261871 phi=0.335357
    pts(12,0)=0.08681635;	 pts(12,1)=0.07016598;	 pts(12,2)=0.81236670;	 //r=0.820000 theta=0.679728 phi=0.136553
    pts(13,0)=0.06802166;	 pts(13,1)=0.05980750;	 pts(13,2)=0.84516041;	 //r=0.850000 theta=0.721227 phi=0.106762
    pts(14,0)=0.32002926;	 pts(14,1)=0.24519047;	 pts(14,2)=0.74830669;	 //r=0.850000 theta=0.653757 phi=0.494174
    return pts;
}
