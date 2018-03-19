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

int cylindre (char namesurf[], char namepatches[], char namepatchcount[], double L, double R, double dt, int*E, int*Nteta, int*Nz, double Ea, double Eb, double Eb2)
{
    int i, j, g, k, save;                       //counters
    double decal_teta;                       // gap along teta.
    double z = -L/2.0f;
    double ez = Ea, iz = Eb, iz2 = Eb2, eo = 1.85f/9.f*2.f*Pi, io = Pi/2.f-eo;    // size of electrode contacts
    int nteta, nl, c, nez, niz, niz2, neo, nio;        //number of points
    double dteta, dl, alpha, o, dez, diz, diz2, deo, dio;
    int np = 0, nt = 0;                //number of points and triangles to create

    decoupe(dt, ez, &nez, &dez);
    decoupe(dt, iz, &niz, &diz);
    decoupe(dt, iz2, &niz2, &diz2);
    decoupe(dt/R, eo, &neo, &deo);
    decoupe(dt/R, io, &nio, &dio);

    nteta = 2*(neo+nio);                //number of points on an outer circle
    dteta = 2*Pi/nteta;
    dt  =  dteta*R;
    decoupe(dt, (L-3*ez-3*iz - iz2)/2, &nl, &dl);
    decal_teta =  -deo*neo/2;
    *Nteta = 4*(neo+nio);
    *Nz = 3*nez +3*niz +niz2+1;

    c = (int)(Pi/2*R/dt+0.5)+1;                   //number of inner circles
    c = (c<2)?2:c;
    alpha = Pi/2/(c-1);                           //distance between inner circles
    std::vector<int> N(c);                      //number of points on inner circles
    std::vector<int> num(c);                    //index of first point on inner circle

    std::vector<double> ang(c);

    int max = (int) floor(nteta*(L/dt+2*c)*4);
    Vertices P(max);
    Triangles T(2*max);

    // indices and coordinates of points before electrode
    for (j=0;j<nl;j++)
    {
        for (i=0;i<nteta;i++)
        {
            o=i*dteta + (j-nl+1)*dteta/2 + decal_teta;
            P[np].x()=R*cos(o);
            P[np].y()=R*sin(o);
            P[np].z()=j*dl +z;
            np ++;
        }
    }

    // indices and coordinates of electrode points
    double zz=nl*dl +z -diz;
    for (j=0;j<=6;j++) {
        int gmax=(j%2 == 0)?((j  ==  4)?niz2:niz):(nez);
        if (j == 0) gmax++;
        for (g=0;g<gmax;g++) {
            zz += (j%2 == 0)?((j  ==  4)?diz2:diz):(dez);
            for (k=0;k<4;k++) {
                for (i=0;i<neo;i++) {
                    o=i*deo +k*(deo*neo+dio*nio) + decal_teta;
                    P[np].x()=R*cos(o);
                    P[np].y()=R*sin(o);
                    P[np].z()=zz;
                    np ++;
                }
                for (i=0;i<nio;i++)
                {
                    o=i*dio +neo*deo +k*(deo*neo+dio*nio) + decal_teta;
                    P[np].x()=R*cos(o);
                    P[np].y()=R*sin(o);
                    P[np].z()=zz;
                    np ++;
                }
            }
        }
    }

    // indices and coordinates of points after electrode
    save=np;
    for (j=nl;j<(2*nl);j++) {
        for (i=0;i<nteta;i++) {
            o=i*dteta + (j-nl)*dteta/2 + decal_teta;
            P[np].x()=R*cos(o);
            P[np].y()=R*sin(o);
            P[np].z()=(j+1)*dl +z +3*ez+3*iz+iz2;
            np ++;
        }
    }


    // definition of triangles before electrode
    for (j=0;j<nl-1;j++) {
        for (i=0;i<nteta;i++) {
            T[nt]=Triangle(P[i+nteta*j], P[(i+1)%nteta +nteta*j], P[i+nteta*(j+1)]);
            T[nt+1]=Triangle(P[i+nteta*j], P[i+nteta*(j+1)], P[(nteta+i-1)%nteta +nteta*(j+1)]);
            nt += 2;
        }
    }

    // definition of transition triangles
    for (i=0;i<nteta;i++) {
        T[nt]=Triangle(P[i+nteta*(nl-1)], P[(2*i+1) + nteta*nl], P[2*i + nteta*nl]);
        nt++;
        T[nt]=Triangle(P[i+nteta*(nl-1)], P[(i+1)%nteta + nteta*(nl-1)], P[2*i+1 + nteta*nl]);
        nt++;
        T[nt]=Triangle(P[(i+1)%nteta+nteta*(nl-1)], P[(2*i+2)%(2*nteta) + nteta*nl], P[(2*i+1) + nteta*nl]);
        nt++;
    }

    // definition of electrode triangles
    for (j=0;j<(3*nez+3*niz+niz2);j++) {
        for (i=0;i<2*nteta;i++) {
            T[nt]=Triangle(P[i+2*nteta*j +nteta*nl], P[i+2*nteta*(j+1) +nteta*nl], P[(2*nteta+i-1)%(2*nteta) +2*nteta*(j+1) +nteta*nl]);
            T[nt+1]=Triangle(P[i+2*nteta*j +nteta*nl], P[(i+1)%(2*nteta)+2*nteta*j +nteta*nl], P[i+2*nteta*(j+1) +nteta*nl]);
            nt += 2;
        }
    }

    // definition of transition triangles
    for (i=0;i<nteta;i++) {
        T[nt]=Triangle(P[i +save], P[2*i  +save-2*nteta], P[(2*i+1) +save-2*nteta]);
        nt++;
        T[nt]=Triangle(P[i +save], P[(2*i+1) +save-2*nteta], P[(i+1)%nteta +save]);
        nt++;
        T[nt]=Triangle(P[(i+1)%nteta +save], P[(2*i+1) +save-2*nteta], P[(2*i+2)%(2*nteta) +save-2*nteta]);
        nt++;
    }

    // definition of triangles after electrode
    for (j=0;j<nl-1;j++) {
        for (i=0;i<nteta;i++) {
            T[nt]=Triangle(P[i+nteta*j +save], P[(i+1)%nteta +nteta*j +save], P[i+nteta*(j+1) +save]);
            T[nt+1]=Triangle(P[i+nteta*j +save], P[i+nteta*(j+1) +save], P[(nteta+i-1)%nteta +nteta*(j+1) +save]);
            nt += 2;
        }
    }

    save = np-nteta;
    // the lids:
    for (g=0;g<=1;g++) {
        num[c-1]=((g == 0)?(0):(save));
        for (j=0;j<c;j++) {
            if (j == 0)
                N[0]=1;
            else
                decoupe(dt/(R*sin(alpha*j)), 2*Pi, &N[j], &ang[j]);
            if (j == c-1) break;
            num[j]=np;
            for (i=0;i<N[j];i++)
            {
                o=i*ang[j] + ((g == 0)?((-nl+1)*dteta/2):((nl-1)*dteta/2)) + decal_teta;
                P[np].x()=R*sin(alpha*j)*cos(o);
                P[np].y()=R*sin(alpha*j)*sin(o);
                P[np].z()=g*L +z + R*cos(alpha*j) *(g*2-1) *0.3;
                np ++;
            }
        }
        for (j=2;j<c;j++) {
            k=0;
            T[nt]=Triangle(P[num[j]], P[num[j]+1], P[num[j-1]]);
            if (g == 0) {
                T[nt].flip();
            }
            nt ++;
            for (i=1;i<N[j];i++)
            {
                if (ang[j-1]*(k+0.5) < ang[j]*(i+0.5))
                {
                    T[nt]=Triangle(P[num[j]+i], P[num[j-1]+(k+1)%N[j-1]], P[num[j-1]+k]);
                    if (g == 0) { 
                        T[nt].flip(); 
                    }
                    k=(k+1)%N[j-1];
                    nt ++;
                }
                T[nt]=Triangle(P[num[j]+i%N[j]], P[num[j]+(i+1)%N[j]], P[num[j-1]+k]);
                if (g == 0) { 
                    T[nt].flip();
                }
                nt ++;
            }
        }
        for (i=0;i<N[1];i++) {
            T[nt]=Triangle(P[num[1]+i], P[num[1]+(i+1)%N[1]], P[num[0]]);
            if (g == 0) {
                T[nt].flip();
            }
            nt ++;
        }
    }
    Mesh surf(P);
    for (i=0;i<nt;i++) {
        surf.push_back(T[i]);
    }
    surf.save(namesurf);

    //the electrode
    if (E[0]>=0){ 
        int trianglecounter=0;
        int electrodetriangleindex = 0;
        int maxelectrodetriangles=0;
        int electrodetriangles [12][20];
        Vect3 electrodecenter;
        for (i=0;i<(nl-1)*2*nteta+3*nteta+1;i++) {
            trianglecounter++;
        }     
        for (j=0;j<=6;j++){
            for (g=0;g<((j%2 == 0)?((j  ==  4)?niz2:niz):(nez));g++) {
                for (k=0;k<4;k++){
                    maxelectrodetriangles = std::max(maxelectrodetriangles, electrodetriangleindex);
                    electrodetriangleindex = 0;
                    for (i=0;i<neo;i++){
                        // find for each electrode to which triangles it corresponds
                        if(j%2 == 1){ 
                            //  first triangle of a pair
                            electrodetriangles[4*(j/2)+k][electrodetriangleindex]=trianglecounter;
                            electrodetriangleindex++;
                            trianglecounter++;
                            //  second triangle of a pair
                            electrodetriangles[4*(j/2)+k][electrodetriangleindex]=trianglecounter;
                            electrodetriangleindex++;
                            trianglecounter++;
                        }
                        else{
                            trianglecounter++;
                            trianglecounter++;
                        }
                    }
                    for (i=0;i<nio;i++)
                    {
                        trianglecounter++;
                        trianglecounter++;
                    }
                }
                //    if(j%2 == 1) eleccounter++;
            }
        }
        for (i= (nl-1)*2*nteta +3*nteta + 2*2*nteta*(3*nez+3*niz+niz2)+1;i<nt;i++){
            trianglecounter++;
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
        for (j=0;j<12;j++){
            for (i=0; i<maxelectrodetriangles;++i){ // compute the barycenters of the triangles corresponding to each electrode
                electrodecenter = (T[electrodetriangles[j][i]].s1() + T[electrodetriangles[j][i]].s2() + T[electrodetriangles[j][i]].s3()) / 3.;
                ofs1 << electrodecenter.x() << ' ' << electrodecenter.y() << ' ' << electrodecenter.z() << std::endl;
            }
            ofs2 << i << std::endl; // the number of patches used to represent electrode j
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
    std::vector<int> Nteta(Nc);
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

        for (i=0;i<Nc;++i) {
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
        for (i=0;i<Nc;++i) {
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
        for (i=0;i<Nc;++i)
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
    for (i=0;i<Nc;i++){
        cylindre(argv[5+i], argv[7], argv[8], L[i], R[i], dt[i], (i == Nc-1)?(Elec):(E), &Nteta[i], &Nz[i], Ea, Eb, Eb2);
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
