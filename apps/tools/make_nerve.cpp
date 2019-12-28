/* FILE: $Id: make_nerve.cpp 2008-04-03 $*/

/*
  Project Name : OpenMEEG

  author            : $Author: mclerc$
  version           : $Revision$
  last revision     : $Date$
  modified by       : $LastChangedBy$
  last modified     : $LastChangedDate$

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

#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cmath>

#include "mesh.h"
#include "sparse_matrix.h"
#include "fast_sparse_matrix.h"
#include "stdlib.h"
#include "triangle.h"
#include "options.h"
#include "geometry.h"
#include "vector.h"

#define Pi 3.14159265358979323846f

using namespace std;
using namespace OpenMEEG;

void getHelp(char** argv);

void decoupe (double dt, double x, int* nx, double* dx)
{
    *nx = (int)(x/dt + 0.5);
    *nx = (*nx > 0)?*nx:1;
    *dx = x/ *nx;
}

int
cylindre(const char namesurf[],const char namepatches[],const char namepatchcount[],
         const double L,const double R,const double d,
         int* E,int* Ntheta,int* Nz,
         const double Ea,const double Eb,const double Eb2)
{
    int  nl, c, nez, niz, niz2, neo, nio;        //number of points
    double dl,  dez, diz, diz2, deo, dio;

    decoupe(d, Ea, &nez, &dez);
    decoupe(d, Eb, &niz, &diz);
    decoupe(d, Eb2, &niz2, &diz2);

    // Size of electrode contacts.

    const double eo = 1.85f/9.f*2.f*Pi;
    const double io = Pi/2.f-eo;

    decoupe(d/R, eo, &neo, &deo);
    decoupe(d/R, io, &nio, &dio);

    const int    ntheta = 2*(neo+nio);                //number of points on an outer circle
    const double dtheta = 2*Pi/ntheta;
    const double dt     =  dtheta*R;
    decoupe(dt, (L-3*Ea-3*Eb - Eb2)/2, &nl, &dl);
    const double decal_theta = -deo*neo/2; // Gap along theta
    *Ntheta = 4*(neo+nio);
    *Nz = 3*nez +3*niz +niz2+1;

    c = (int)(Pi/2*R/dt+0.5)+1;                   //number of inner circles
    c = (c<2)?2:c;
    const double alpha = Pi/2/(c-1);                           //distance between inner circles
    std::vector<int> N(c);                      //number of points on inner circles
    std::vector<int> num(c);                    //index of first point on inner circle

    std::vector<double> ang(c);

    int max = (int) floor(ntheta*(L/dt+2*c)*4);
    Vertices P(max);
    Triangles T(2*max);

    // indices and coordinates of points before electrode

    unsigned np = 0;    //  Number of points.
    const double z = -L/2.0f;
    for (unsigned j=0;j<nl;j++) {
        for (unsigned i=0;i<ntheta;i++) {
            const double angle = i*dtheta+(j-nl+1)*dtheta/2+decal_theta;
            P[np].x()=R*cos(angle);
            P[np].y()=R*sin(angle);
            P[np].z()=j*dl +z;
            np ++;
        }
    }

    // indices and coordinates of electrode points
    double zz=nl*dl +z -diz;
    for (unsigned j=0;j<=6;j++) {
        int gmax=(j%2 == 0)?((j  ==  4)?niz2:niz):(nez);
        if (j == 0) gmax++;
        for (unsigned g=0;g<gmax;g++) {
            zz += (j%2 == 0)?((j  ==  4)?diz2:diz):(dez);
            for (unsigned k=0;k<4;k++) {
                for (unsigned i=0;i<neo;i++) {
                    const double angle=i*deo +k*(deo*neo+dio*nio) + decal_theta;
                    P[np].x()=R*cos(angle);
                    P[np].y()=R*sin(angle);
                    P[np].z()=zz;
                    np ++;
                }
                for (unsigned i=0;i<nio;i++)
                {
                    const double angle=i*dio +neo*deo +k*(deo*neo+dio*nio) + decal_theta;
                    P[np].x()=R*cos(angle);
                    P[np].y()=R*sin(angle);
                    P[np].z()=zz;
                    np ++;
                }
            }
        }
    }

    // indices and coordinates of points after electrode
    unsigned save=np;
    for (unsigned j=nl;j<(2*nl);j++) {
        for (unsigned i=0;i<ntheta;i++) {
            const double angle = i*dtheta+(j-nl)*dtheta/2+decal_theta;
            P[np].x() = R*cos(angle);
            P[np].y() = R*sin(angle);
            P[np].z() = (j+1)*dl+z+3*Ea+3*Eb+Eb2;
            ++np;
        }
    }

    // definition of triangles before electrode

    unsigned nt = 0;  // Number of triangles.
    for (unsigned j=0;j<nl-1;j++)
        for (unsigned i=0;i<ntheta;i++) {
            T[nt++] = Triangle(P[i+ntheta*j], P[(i+1)%ntheta +ntheta*j], P[i+ntheta*(j+1)]);
            T[nt++] = Triangle(P[i+ntheta*j], P[i+ntheta*(j+1)], P[(ntheta+i-1)%ntheta +ntheta*(j+1)]);
        }

    // definition of transition triangles
    for (unsigned i=0;i<ntheta;i++) {
        T[nt++] = Triangle(P[i+ntheta*(nl-1)], P[(2*i+1) + ntheta*nl], P[2*i + ntheta*nl]);
        T[nt++] = Triangle(P[i+ntheta*(nl-1)], P[(i+1)%ntheta + ntheta*(nl-1)], P[2*i+1 + ntheta*nl]);
        T[nt++] = Triangle(P[(i+1)%ntheta+ntheta*(nl-1)], P[(2*i+2)%(2*ntheta) + ntheta*nl], P[(2*i+1) + ntheta*nl]);
    }

    // definition of electrode triangles
    for (unsigned j=0;j<(3*nez+3*niz+niz2);j++)
        for (unsigned i=0;i<2*ntheta;i++) {
            T[nt++] = Triangle(P[i+2*ntheta*j+ntheta*nl],P[i+2*ntheta*(j+1)+ntheta*nl],P[(2*ntheta+i-1)%(2*ntheta)+2*ntheta*(j+1)+ntheta*nl]);
            T[nt++] = Triangle(P[i+2*ntheta*j+ntheta*nl],P[(i+1)%(2*ntheta)+2*ntheta*j+ntheta*nl],P[i+2*ntheta*(j+1)+ntheta*nl]);
        }

    // definition of transition triangles
    for (unsigned i=0;i<ntheta;i++) {
        T[nt++] = Triangle(P[i+save],P[2*i+save-2*ntheta],P[(2*i+1)+save-2*ntheta]);
        T[nt++] = Triangle(P[i+save],P[(2*i+1)+save-2*ntheta],P[(i+1)%ntheta+save]);
        T[nt++] = Triangle(P[(i+1)%ntheta+save],P[(2*i+1)+save-2*ntheta],P[(2*i+2)%(2*ntheta)+save-2*ntheta]);
    }

    // definition of triangles after electrode
    for (unsigned j=0;j<nl-1;j++) {
        for (unsigned i=0;i<ntheta;i++) {
            T[nt++] = Triangle(P[i+ntheta*j+save],P[(i+1)%ntheta+ntheta*j+save],P[i+ntheta*(j+1)+save]);
            T[nt++] = Triangle(P[i+ntheta*j+save],P[i+ntheta*(j+1)+save],P[(ntheta+i-1)%ntheta+ntheta*(j+1)+save]);
        }
    }

    save = np-ntheta;
    // the lids:
    for (unsigned g=0;g<=1;g++) {
        num[c-1]=((g == 0)?(0):(save));
        for (unsigned j=0;j<c;j++) {
            if (j==0)
                N[0]=1;
            else
                decoupe(dt/(R*sin(alpha*j)), 2*Pi, &N[j], &ang[j]);
            if (j==c-1) break;
            num[j] = np;
            for (unsigned i=0;i<N[j];++i,++np) {
                const double angle = i*ang[j]+((g==0) ? 1-nl : nl-1)*dtheta/2+decal_theta;
                P[np].x() = R*sin(alpha*j)*cos(angle);
                P[np].y() = R*sin(alpha*j)*sin(angle);
                P[np].z() = g*L+z+R*cos(alpha*j)*(g*2-1)*0.3;
            }
        }
        for (unsigned j=2;j<c;j++) {
            unsigned k=0;
            T[nt] = Triangle(P[num[j]],P[num[j]+1],P[num[j-1]]);
            if (g==0)
                T[nt].change_orientation();
            ++nt;
            for (unsigned i=1;i<N[j];i++) {
                if (ang[j-1]*(k+0.5)<ang[j]*(i+0.5)) {
                    T[nt] = Triangle(P[num[j]+i], P[num[j-1]+(k+1)%N[j-1]], P[num[j-1]+k]);
                    if (g==0)
                        T[nt].change_orientation();
                    k = (k+1)%N[j-1];
                    ++nt;
                }
                T[nt] = Triangle(P[num[j]+i%N[j]], P[num[j]+(i+1)%N[j]], P[num[j-1]+k]);
                if (g==0)
                    T[nt].change_orientation();
                ++nt;
            }
        }
        for (unsigned i=0;i<N[1];++i,++nt) {
            T[nt] = Triangle(P[num[1]+i],P[num[1]+(i+1)%N[1]],P[num[0]]);
            if (g==0)
                T[nt].change_orientation();
        }
    }
    Mesh surf(std::move(P),std::move(T));
    surf.save(namesurf);

    //the electrode
    if (E[0]>=0){ 
        unsigned maxelectrodetriangles=0;
        unsigned electrodetriangles[12][20];

        for (unsigned j=0,electrodetriangleindex=0,trianglecounter=(nl-1)*2*ntheta+3*ntheta+1; j<=6; ++j) {
            const unsigned ng = (j%2==0) ? ((j==4) ? niz2 : niz) : nez;
            for (unsigned g=0; g<ng; ++g) {
                for (unsigned k=0; k<4; ++k,trianglecounter+=2*nio) {
                    maxelectrodetriangles = std::max(maxelectrodetriangles, electrodetriangleindex);
                    for (unsigned i=0; i<neo; ++i,trianglecounter+=2){
                        // find for each electrode to which triangles it corresponds
                        if (j%2==1) { 
                            electrodetriangles[4*(j/2)+k][electrodetriangleindex++] = trianglecounter;
                            electrodetriangles[4*(j/2)+k][electrodetriangleindex++] = trianglecounter+1;
                        }
                    }
                }
                //    if(j%2 == 1) eleccounter++;
            }
        }

        // compute and store the electrodes, represented as a collection of 3D positions (patches)

        // A file for storing the patch coordinates (for all the electrodes)

        std::ofstream ofs1(namepatches);
        if (!ofs1) {
            std::cerr << "Cannot open file " << namepatches << " for writing." << std::endl;
            exit(1);
        }

        // A file for storing the number of patches per electrode (12 lines)
        std::ofstream ofs2(namepatchcount);
        if (!ofs2) {
            std::cerr << "Cannot open file " << namepatchcount << " for writing." << std::endl;
            exit(1);
        }
        Vect3 electrodecenter;
        for (unsigned j=0;j<12;j++) {
            for (unsigned i=0; i<maxelectrodetriangles;++i){ // compute the barycenters of the triangles corresponding to each electrode
                electrodecenter = (T[electrodetriangles[j][i]].vertex(0)+T[electrodetriangles[j][i]].vertex(1)+T[electrodetriangles[j][i]].vertex(2))/3;
                ofs1 << electrodecenter.x() << ' ' << electrodecenter.y() << ' ' << electrodecenter.z() << std::endl;
            }
            ofs2 << maxelectrodetriangles << std::endl; // the number of patches used to represent electrode j
        }
    }
    return 0;
}

template <typename T>
void input(const char* message,T& value) {
    std::cout << message;
    std::cout.flush();
    std::cin >> value;
}

int main(int argc, char** argv)
{
    print_version(argv[0]);

    command_usage("Make nerve geometry from existing parameters or make nerve geometry and parameter file from commandline user interface.");
    if ((!strcmp(argv[1], "-h")) | (!strcmp(argv[1], "--help"))) getHelp(argv);
    disp_argv(argc, argv);
    // char nom[80], s[80];
    const int Nc = 2;
    int i, Elec[4], E[1]={-1};
    double Ea, Eb, Eb2;
    std::vector<int> Ntheta(Nc);
    std::vector<int> Nz(Nc);
    std::vector<double> L(Nc);
    std::vector<double> R(Nc);
    std::vector<double> dt(Nc);
    std::vector<double> sig(Nc+1);
    sig[Nc]=0;

    if (!strcmp(argv[1], "-makeparameters")){
        // input("number of cylinders (for the moment fixed to 2) :",Nc);
        input("Size of electrode contact along z :",Ea);
        input("Distance between first anode and cathode:",Eb);
        input("Distance between second anode and cathode :",Eb2);

        for (unsigned i=0;i<Nc;++i) {
            std::cout << " Nerve geometry (" << ((i==0) ? "inner cylinder" : "outer cylinder containing electrode") << std::endl;
            input("Length :",L[1]);
            input("Radius :",R[i]);
            input("Discretization step :",dt[i]);
            input("Conductivity :",sig[i]);
        }

        std::ofstream ofs(argv[2]);
        if (!ofs) {
            std::cerr << "Cannot open file " << argv[2] << " for writing." << std::endl;
            return 1;
        }
        ofs << Ea << std::endl
            << Eb << std::endl
            << Eb2 << std::endl;
        for (unsigned i=0;i<Nc;++i) {
            ofs << L[i] << std::endl
                << R[i] << std::endl
                << dt[i] << std::endl
                << sig[i] << std::endl;
        }
    }
    else if (!strcmp(argv[1], "-useparameters")){
        // Read parameter file :
        std::ifstream ifs(argv[2]);
        if (!ifs) {
            std::cerr << "Cannot open file " << argv[2] << " for reading." << std::endl;
            return 1;
        }

        ifs >> Ea >> Eb >> Eb2;
        for (unsigned i=0;i<Nc;++i)
            ifs >> L[i] >> R[i] >> dt[i] >> sig[i];
    } else {
        std::cerr << "unknown argument: " << argv[1] << endl;
        return 1;
    }

    //  the electrode:
    Elec[0]=(int)((L[Nc-1]/2-Ea/2-Eb/2)/dt[Nc-1] +0.5);
    Elec[1]=(int)(Eb / dt[Nc-1]);
    Elec[2]=(int)((L[Nc-1]/2+Ea/2-Eb2/2)/dt[Nc-1] +0.5);
    Elec[3]=(int)(Eb2 /dt[Nc-1] +0.5);
    // saving the corresponding geom and cond files

    std::ofstream ofgeom(argv[3]);
    if (!ofgeom) {
        std::cerr << "Cannot open file " << argv[3] << " for writing." << std::endl;
        return 1;
    }

    ofgeom << "# Domain Description 1.0" << std::endl << std::endl
           << "Interfaces " << Nc << " Mesh" << std::endl << std::endl;
    for (unsigned i=0;i<Nc;i++){
        cylindre(argv[5+i], argv[7], argv[8], L[i], R[i], dt[i], (i == Nc-1)?(Elec):(E), &Ntheta[i], &Nz[i], Ea, Eb, Eb2);
        const std::string& name = argv[5+i];
        // only keep the file name without the directories (after the last /)
        ofgeom << name.substr(name.rfind("/")+1) << std::endl;
    }
    ofgeom << std::endl
           << "Domains " << Nc+1 << std::endl << std::endl
           << "Domain CSF 1 -2" << std::endl
           << "Domain Nerve -1 shared" << std::endl
           << "Domain Air 2" << std::endl;

    std::ofstream ofcond(argv[4]);
    if (!ofcond) {
        std::cerr << "Cannot open file " << argv[4] << " for writing." << std::endl;
        return 2;
    }
    ofcond << "# Properties Description 1.0 (Conductivities)" << std::endl << std::endl
           << "Air        0.0" << std::endl
           << "CSF        " << sig[1] << std::endl
           << "Nerve      " << sig[0] << std::endl;

    std::cerr << "New nerve geometry made." << std::endl;

    return 0;
}

void getHelp(char** argv) {
    std::cerr << argv[0] <<" [options] [filepaths...]"     << std::endl << std::endl
              << "       Arguments :"                      << std::endl << std::endl
              << "         -useparameters       "          << std::endl
              << "               Input parameter file"     << std::endl
              << "               Output geom file "        << std::endl
              << "               Output cond file "        << std::endl
              << "               Output surf1 tri mesh "   << std::endl
              << "               Output surf2 tri mesh "   << std::endl
              << "               Output patches Matrix"    << std::endl
              << "               Output patchcount Matrix" << std::endl << std::endl
              << "         -makeparameters       "         << std::endl
              << "               Output parameter file"    << std::endl
              << "               Output geom file "        << std::endl
              << "               Output cond file "        << std::endl
              << "               Output surf1 tri mesh "   << std::endl
              << "               Output surf2 tri mesh "   << std::endl
              << "               Output patches Matrix"    << std::endl
              << "               Output patchcount Matrix" << std::endl << std::endl;
    exit(0);
}
