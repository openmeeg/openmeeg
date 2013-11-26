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

#include <fstream>
#include <cstring>
#include <vector>

#include "mesh.h"
#include "sparse_matrix.h"
#include "fast_sparse_matrix.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "triangle.h"
#include "options.h"
#include "geometry.h"
#include "vector.h"

#define Pi 3.14159265358979323846f

using namespace std;
using namespace OpenMEEG;

void getHelp(char** argv);

void decoupe (double dt, double x, int*nx, double*dx)
{
    *nx = (int)(x/dt + 0.5);
    *nx = (nx > 0)?*nx:1;
    *dx = x/ *nx;
}

int cylindre (char namesurf[], char namepatches[], char namepatchcount[], double L, double R, double dt, int*E, int*Nteta, int*Nz, double Ea, double Eb, double Eb2)
{
    FILE *F, *G;
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
        F = fopen(namepatches, "w");   // a file for storing the patch coordinates (for all the electrodes)
        G = fopen(namepatchcount, "w"); // a file for storing the number of patches per electrode (12 lines)
        for (j=0;j<12;j++){
            for (i=0; i<maxelectrodetriangles;i++){ // compute the barycenters of the triangles corresponding to each electrode
                electrodecenter = (T[electrodetriangles[j][i]].s1() + T[electrodetriangles[j][i]].s2() + T[electrodetriangles[j][i]].s3()) / 3.;
                fprintf(F, "%f %f %f\n", electrodecenter.x(), electrodecenter.y(), electrodecenter.z());	
            }
            fprintf(G, "%d\n", i); // the number of patches used to represent electrode j
        }
        fclose(F);
        fclose(G);
    }
    return 0;
}
int main(int argc, char** argv)
{
    print_version(argv[0]);

    command_usage("Make nerve geometry from existing parameters or make nerve geometry and parameter file from commandline user interface.");
    if ((!strcmp(argv[1], "-h")) | (!strcmp(argv[1], "--help"))) getHelp(argv);
    disp_argv(argc, argv);
    FILE *F,*Fgeom,*Fcond;
    // char nom[80], s[80];
    int i, Nc=2, Elec[4], E[1]={-1};
    double Ea, Eb, Eb2;
    std::vector<int> Nteta(Nc);
    std::vector<int> Nz(Nc);
    std::vector<double> L(Nc);
    std::vector<double> R(Nc);
    std::vector<double> dt(Nc);
    std::vector<double> sig(Nc+1);
    string name;
    sig[Nc]=0;

    if (!strcmp(argv[1], "-makeparameters")){
        // fprintf(stdout, "number of cylinders (for the moment fixed to 2) :");fflush(stdout);
        // fscanf(stdin, "%d", &Nc);
        fprintf(stdout, "size of electrode contact along z :");fflush(stdout);
        fscanf(stdin, "%lf", &Ea);
        fprintf(stdout, "distance between first anode and cathode:");fflush(stdout);
        fscanf(stdin, "%lf", &Eb);
        fprintf(stdout, "distance between second anode and cathode :");fflush(stdout);
        fscanf(stdin, "%lf", &Eb2);

        for (i=0;i<Nc;i++)
        {
            if (i == 0)
            {printf(" nerve geometry (inner cylinder)\n");}
            if (i == 1)
            {printf(" nerve geometry (outer cylinder containing electrode)\n");}
            fprintf(stdout, "length :");fflush(stdout);
            fscanf(stdin, "%lf", &L[i]);
            fprintf(stdout, "radius :");fflush(stdout);
            fscanf(stdin, "%lf", &R[i]);
            fprintf(stdout, "discretization step :");fflush(stdout);
            fscanf(stdin, "%lf", &dt[i]);
            fprintf(stdout, "conductivity :");fflush(stdout);
            fscanf(stdin, "%lf", &sig[i]);
        }
        F = fopen (argv[2], "w");
        fprintf(F, "%lf\n", Ea);
        fprintf(F, "%lf\n", Eb);
        fprintf(F, "%lf\n", Eb2);
        for (i=0;i<Nc;i++)
        {
            fprintf(F, "%lf\n", L[i]);
            fprintf(F, "%lf\n", R[i]);
            fprintf(F, "%lf\n", dt[i]);
            fprintf(F, "%lf\n", sig[i]);
        }
        fclose(F);
    }
    else if (!strcmp(argv[1], "-useparameters")){
        // Read parameter file :
        F = fopen (argv[2], "r");
        if (F == NULL) {perror ("error fopen: \n"); return -1;}
        fscanf(F, "%lf", &Ea);
        fscanf(F, "%lf", &Eb);;
        fscanf(F, "%lf", &Eb2);
        for (i=0;i<Nc;i++)
        {
            fscanf(F, "%lf", &L[i]);
            fscanf(F, "%lf", &R[i]);
            fscanf(F, "%lf", &dt[i]);
            fscanf(F, "%lf", &sig[i]);
            /*      fclose(F);*/
        }
        fclose(F);
    }

    else  cerr << "unknown argument: " << argv[1] << endl;
    //  the electrode:
    Elec[0]=(int)((L[Nc-1]/2-Ea/2-Eb/2)/dt[Nc-1] +0.5);
    Elec[1]=(int)(Eb / dt[Nc-1]);
    Elec[2]=(int)((L[Nc-1]/2+Ea/2-Eb2/2)/dt[Nc-1] +0.5);
    Elec[3]=(int)(Eb2 /dt[Nc-1] +0.5);
    // saving the corresponding geom and cond files
    Fgeom = fopen (argv[3], "w");
    Fcond = fopen (argv[4], "w");
    if (Fgeom == NULL) {perror ("error fopen: .geom\n"); return -1;}
    if (Fcond == NULL) {perror ("error fopen: .cond\n"); return -1;}
    fprintf(Fgeom, "# Domain Description 1.0\n\n");
    fprintf(Fgeom, "Interfaces %d Mesh\n\n", Nc);
    fprintf(Fcond, "# Properties Description 1.0 (Conductivities)\n\n");
    for (i=0;i<Nc;i++){
        cylindre(argv[5+i], argv[7], argv[8], L[i], R[i], dt[i], (i == Nc-1)?(Elec):(E), &Nteta[i], &Nz[i], Ea, Eb, Eb2);
        name = argv[5+i];
        name = name.substr(name.rfind("/")+1); // only keep the file name without the directories (after the last /)
        fprintf(Fgeom, "%s\n", name.c_str());    // c_str to convert string back to char
    }
    fprintf(Fcond, "Air        0.0\n");
    fprintf(Fcond, "CSF        %lf\n", sig[1]);
    fprintf(Fcond, "Nerve      %lf\n", sig[0]);
    fprintf(Fgeom, "\nDomains %d\n\n", Nc+1);
    fprintf(Fgeom, "Domain CSF 1 -2\n");
    fprintf(Fgeom, "Domain Nerve -1 shared\n");
    fprintf(Fgeom, "Domain Air 2\n");
    fclose(Fgeom);
    fclose(Fcond);
    printf("New nerve geometry made.\n");
    return 0;
}

void getHelp(char** argv) {
    cout << argv[0] <<" [options] [filepaths...]" << endl << endl;
    cout << "       Arguments :" << endl <<endl;
    cout << "         -useparameters       " <<endl;
    cout << "               Input parameter file" << endl;
    cout << "               Output geom file " << endl;
    cout << "               Output cond file " << endl;
    cout << "               Output surf1 tri mesh " << endl;
    cout << "               Output surf2 tri mesh " << endl;
    cout << "               Output patches Matrix" << endl;
    cout << "               Output patchcount Matrix" << endl << endl;
    cout << "         -makeparameters       " <<endl;
    cout << "               Output parameter file" << endl;
    cout << "               Output geom file " << endl;
    cout << "               Output cond file " << endl;
    cout << "               Output surf1 tri mesh " << endl;
    cout << "               Output surf2 tri mesh " << endl;
    cout << "               Output patches Matrix" << endl;
    cout << "               Output patchcount Matrix" << endl << endl;
    exit(0);
}
