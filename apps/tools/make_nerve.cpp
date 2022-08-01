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
#include <functional>

#include "mesh.h"
#include "sparse_matrix.h"
#include "fast_sparse_matrix.h"
#include "stdlib.h"
#include "triangle.h"
#include "commandline.h"
#include "constants.h"
#include "geometry.h"
#include "vector.h"

using namespace OpenMEEG;

void getHelp(char** argv);

unsigned num_points(const double x,const double dt) { return std::max(static_cast<unsigned>(x/dt+0.5),1U); }
double   step_size(const double x,const unsigned n) { return x/n; }

typedef std::vector<TriangleIndices> MeshTriangles;

template <typename T>
void input(const char* message,T& value) {
    std::cout << message;
    std::cout.flush();
    std::cin >> value;
}

void make_lid(const unsigned n,const double R,const double z0,const double dt,const double theta0,Vertices& vertices,MeshTriangles& triangles) {

    const unsigned c     = std::max(2U,static_cast<unsigned>((Pi/2*R/dt+0.5)+1)); // Number of inner circles
    const double   alpha = Pi/2/(c-1);                                            // Distance between inner circles

    std::vector<unsigned> N(c);   // Number of points on inner circles
    std::vector<unsigned> num(c); // Index of first point on inner circle

    std::vector<double> ang(c);

    N[0] = 1;
    ang[0] = 2*Pi;
    for (unsigned j=1; j<c; ++j) {
        N[j]   = num_points(2*Pi,dt/(R*sin(alpha*j)));
        ang[j] = step_size(2*Pi,N[j]);
    }

    const double dz = (n==0) ? -0.3 : 0.3;

    for (unsigned j=0; j<c-1; ++j) {
        num[j] = vertices.size();
        for (double angle=theta0; angle<2*Pi; angle+=ang[j])
            vertices.push_back(Vertex(R*sin(alpha*j)*cos(angle),R*sin(alpha*j)*sin(angle),z0+R*cos(alpha*j)*dz));
    }
    num[c-1] = n;

    std::function<void(const unsigned i,const unsigned j,const unsigned k)> add_triangle;
    if (n==0) {
        add_triangle = [&](const unsigned i,const unsigned j,const unsigned k) { triangles.push_back({ i, j, k }); };
    } else {
        add_triangle = [&](const unsigned i,const unsigned j,const unsigned k) { triangles.push_back({ j, i, k }); };
    }

    for (unsigned j=2; j<c; ++j) {
        unsigned k=0;
        add_triangle(num[j],num[j]+1,num[j-1]);
        for (unsigned i=1; i<N[j]; ++i) {
            if (ang[j-1]*(k+0.5)<ang[j]*(i+0.5)) {
                add_triangle(num[j]+i,num[j-1]+(k+1)%N[j-1],num[j-1]+k);
                k = (k+1)%N[j-1];
            }
            add_triangle(num[j]+i%N[j],num[j]+(i+1)%N[j],num[j-1]+k);
        }
    }

    for (unsigned i=0; i<N[1]; ++i)
        add_triangle(num[1]+i,num[1]+(i+1)%N[1],num[0]);
}

void
cylindre(const std::string& namesurf,const char namepatches[],const char namepatchcount[],
         const double L,const double R,const double d,const unsigned* E,
         unsigned& Ntheta,unsigned& Nz,
         const double Ea,const double Eb,const double Eb2)
{
    const unsigned nez = num_points(Ea,d);
    const double   dez = step_size(Ea,nez);
    const unsigned niz = num_points(Eb,d);
    const double   diz = step_size(Eb,niz);
    const unsigned niz2 = num_points(Eb2,d);
    const double   diz2 = step_size(Eb2,niz2);

    // Size of electrode contacts.

    const double eo = 1.85f/9.f*2.f*Pi;
    const double io = Pi/2.f-eo;

    const unsigned neo = num_points(eo,d/R);
    const double   deo = step_size(eo,neo);

    const unsigned nio = num_points(io,d/R);
    const double   dio = step_size(io,nio);

    const unsigned ntheta = 2*(neo+nio);                //number of points on an outer circle
    const double   dtheta = 2*Pi/ntheta;
    const double   dt     =  dtheta*R;

    const double   El = (L-3*Ea-3*Eb-Eb2)/2;
    const unsigned nl = num_points(El,dt);
    const double   dl = step_size(El,nl);

    const double decal_theta = -deo*neo/2; // Gap along theta
    Ntheta = 4*(neo+nio);
    Nz     = 3*nez+3*niz+niz2+1;

    // Points before electrode

    Vertices vertices;

    const double z = -L/2.0f;
    for (unsigned j=0; j<nl; ++j) {
        for (unsigned i=0; i<ntheta; ++i) {
            const double angle = i*dtheta+(j-nl+1)*dtheta/2+decal_theta;
            vertices.push_back(Vertex(R*cos(angle),R*sin(angle),j*dl+z));
        }
    }

    // indices and coordinates of electrode points

    double zz = nl*dl+z-diz;
    const unsigned gmax[] = { 1+niz, nez, niz, nez, niz2, nez, niz };
    const double   dzz[]  = {   diz, dez, diz, dez, diz2, dez, diz };
    for (unsigned j=0,l=1; j<=6; ++j,l=0) {
        for (unsigned g=0; g<gmax[j]; ++g) {
            zz += dzz[j];
            for (unsigned k=0; k<4; ++k) {
                for (unsigned i=0; i<neo; ++i) {
                    const double angle = i*deo+k*(deo*neo+dio*nio)+decal_theta;
                    vertices.push_back(Vertex(R*cos(angle),R*sin(angle),zz));
                }
                for (unsigned i=0; i<nio; ++i) {
                    const double angle = i*dio+neo*deo+k*(deo*neo+dio*nio)+decal_theta;
                    vertices.push_back(Vertex(R*cos(angle),R*sin(angle),zz));
                }
            }
        }
    }

    // indices and coordinates of points after electrode

    const unsigned n1 = vertices.size();
    for (unsigned j=nl; j<(2*nl); ++j) {
        for (unsigned i=0; i<ntheta; ++i) {
            const double angle = i*dtheta+(j-nl)*dtheta/2+decal_theta;
            vertices.push_back(Vertex(R*cos(angle),R*sin(angle),(j+1)*dl+z+3*Ea+3*Eb+Eb2));
        }
    }

    // Triangles before electrode

    MeshTriangles triangles;

    for (unsigned j=0; j<nl-1; ++j)
        for (unsigned i=0; i<ntheta; ++i) {
            const unsigned ind1 = i+ntheta*j;
            const unsigned ind2 = (i+1)%ntheta+ntheta*j;
            const unsigned ind3 = ind1+ntheta;
            const unsigned ind4 = (ntheta+i-1)%ntheta+ind3-i;
            triangles.push_back({ ind1, ind2, ind3 });
            triangles.push_back({ ind1, ind3, ind4 });
        }

    // Transition triangles

    for (unsigned i=0; i<ntheta; ++i) {
        const unsigned ind1 = i+ntheta*(nl-1);
        const unsigned ind3 = 2*i+ntheta*nl;
        const unsigned ind2 = ind3+1;
        const unsigned ind4 = (i+1)%ntheta+ntheta*(nl-1);
        const unsigned ind5 = (2*i+2)%(2*ntheta)+ntheta*nl;
        triangles.push_back({ ind1, ind2, ind3 });
        triangles.push_back({ ind1, ind4, ind2 });
        triangles.push_back({ ind4, ind5, ind2 });
    }

    // Electrode triangles

    for (unsigned j=0; j<(3*nez+3*niz+niz2); ++j)
        for (unsigned i=0; i<2*ntheta; ++i) {
            const unsigned ind1 = i+2*ntheta*j+ntheta*nl;
            const unsigned ind2 = ind1+2*ntheta;
            const unsigned ind3 = (2*ntheta+i-1)%(2*ntheta)+ind2-i;
            const unsigned ind4 = (i+1)%(2*ntheta)+ind1-i;
            triangles.push_back({ ind1, ind2, ind3 });
            triangles.push_back({ ind1, ind4, ind2 });
        }

    // Transition triangles

    for (unsigned i=0; i<ntheta; ++i) {
        const unsigned ind1 = i+n1;
        const unsigned ind2 = ind1+i-2*ntheta;
        const unsigned ind3 = ind2+1;
        const unsigned ind4 = (i+1)%ntheta+n1;
        const unsigned ind5 = (2*i+2)%(2*ntheta)+n1-2*ntheta;
        triangles.push_back({ ind1, ind2, ind3 });
        triangles.push_back({ ind1, ind3, ind4 });
        triangles.push_back({ ind4, ind3, ind5 });
    }

    // Triangles after electrode

    for (unsigned j=0; j<nl-1; ++j) {
        for (unsigned i=0; i<ntheta; ++i) {
            const unsigned ind1 = i+ntheta*j+n1;
            const unsigned ind2 = (i+1)%ntheta+ntheta*j+n1;
            const unsigned ind3 = ind1+ntheta;
            const unsigned ind4 = (ntheta+i-1)%ntheta+ind3-i;
            triangles.push_back({ ind1, ind2, ind3 });
            triangles.push_back({ ind1, ind3, ind4 });
        }
    }

    const unsigned n2 = vertices.size()-ntheta;

    //  Lids:

    make_lid(0,R,z,  dt,(1-nl)*dtheta/2+decal_theta,vertices,triangles);
    make_lid(1,R,z+L,dt,(nl-1)*dtheta/2+decal_theta,vertices,triangles);

    //  Since we created the points, there is no need to map triangle indices.

    Mesh surf;

    for (const auto& vertex : vertices)
        surf.geometry().vertices().push_back(vertex);

    for (const auto& triangle : triangles)
        surf.add_triangle(triangle);

    surf.update(true);
        
    surf.save(namesurf);

    // Electrode

    if (E!=nullptr) { 

        unsigned maxelectrodetriangles = 0;
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
        for (unsigned j=0; j<12; ++j) {
            for (unsigned i=0; i<maxelectrodetriangles; ++i){ // compute the barycenters of the triangles corresponding to each electrode
                electrodecenter = (surf.triangles()[electrodetriangles[j][i]].vertex(0)+surf.triangles()[electrodetriangles[j][i]].vertex(1)+surf.triangles()[electrodetriangles[j][i]].vertex(2))/3;
                ofs1 << electrodecenter.x() << ' ' << electrodecenter.y() << ' ' << electrodecenter.z() << std::endl;
            }
            ofs2 << maxelectrodetriangles << std::endl; // the number of patches used to represent electrode j
        }
    }
}

int main(int argc, char** argv) {

    const CommandLine cmd(argc,argv,"Make nerve geometry from existing parameters or make nerve geometry and parameter file from commandline user interface.");
    if ((!strcmp(argv[1],"-h")) || (!strcmp(argv[1], "--help")))
        getHelp(argv);
    print_version(argv[0]);
    print_commandline(argc,argv);

    // input("number of cylinders (for the moment fixed to 2) :",Nc);

    const unsigned Nc = 2;

    double Ea, Eb, Eb2;
    std::vector<double> L(Nc);
    std::vector<double> R(Nc);
    std::vector<double> dt(Nc);
    std::vector<double> sig(Nc+1);
    sig[Nc]=0;

    if (!strcmp(argv[1], "-makeparameters")){
        input("Size of electrode contact along z :",Ea);
        input("Distance between first anode and cathode:",Eb);
        input("Distance between second anode and cathode :",Eb2);

        for (unsigned i=0; i<Nc; ++i) {
            std::cout << " Nerve geometry (" << ((i==0) ? "inner cylinder" : "outer cylinder containing electrode") << std::endl;
            input("Length :",L[i]);
            input("Radius :",R[i]);
            input("Discretization step :",dt[i]);
            input("Conductivity :",sig[i]);
        }

        std::ofstream ofs(argv[2]);
        if (!ofs) {
            std::cerr << "Cannot open file " << argv[2] << " for writing." << std::endl;
            return 1;
        }
        ofs << Ea << std::endl << Eb << std::endl << Eb2 << std::endl;
        for (unsigned i=0; i<Nc; ++i)
            ofs << L[i] << ' ' << R[i] << ' ' << dt[i] << ' ' << sig[i] << std::endl;

    } else if (!strcmp(argv[1], "-useparameters")){

        // Read parameter file :

        std::ifstream ifs(argv[2]);
        if (!ifs) {
            std::cerr << "Cannot open file " << argv[2] << " for reading." << std::endl;
            return 1;
        }

        ifs >> Ea >> Eb >> Eb2;
        for (unsigned i=0; i<Nc; ++i)
            ifs >> L[i] >> R[i] >> dt[i] >> sig[i];

        std::cerr << L[0] << ' ' << L[1] << std::endl;
    } else {
        std::cerr << "unknown argument: " << argv[1] << std::endl;
        return 1;
    }

    //  the electrode:

    const unsigned Elec[4] = {
        static_cast<unsigned>((L[Nc-1]/2-Ea/2-Eb/2)/dt[Nc-1]+0.5),  static_cast<unsigned>(Eb/dt[Nc-1]),
        static_cast<unsigned>((L[Nc-1]/2+Ea/2-Eb2/2)/dt[Nc-1]+0.5), static_cast<unsigned>(Eb2/dt[Nc-1]+0.5)
    };
    
    // saving the corresponding geom and cond files

    std::ofstream ofgeom(argv[3]);
    if (!ofgeom) {
        std::cerr << "Cannot open file " << argv[3] << " for writing." << std::endl;
        return 1;
    }

    ofgeom << "# Domain Description 1.0" << std::endl << std::endl
           << "Interfaces " << Nc << " Mesh" << std::endl << std::endl;

    std::vector<unsigned> Ntheta(Nc);
    std::vector<unsigned> Nz(Nc);

    for (unsigned i=0;i<Nc;i++) {
        const unsigned* electrode = (i==Nc-1) ? Elec : nullptr;
        const std::string& name = argv[5+i];
        cylindre(name,argv[7],argv[8],L[i],R[i],dt[i],electrode,Ntheta[i],Nz[i],Ea,Eb,Eb2);
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
