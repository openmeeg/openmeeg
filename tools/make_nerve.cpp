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

#include "mesh3.h"
#include "sparse_matrice.h"
#include "fast_sparse_matrice.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "triangle.h"
#include "options.h"
#include "geometry.h"
#include <vecteur.h>
#define Pi acos(-1)
using namespace std;

void getHelp(char** argv);

void decoupe (float dt, float x, int*nx, float*dx)
{
  *nx=(int)(x/dt+0.5);
  *nx= (nx>0)?*nx:1;
  *dx=x/ *nx;
}

void permute (Triangle *t)
{
  int temp= (*t)[1];
  (*t)[1]=(*t)[2];
  (*t)[2]=temp;
}

int cylindre (char namesurf[],char namestim[],float L,float R,float dt,int*E,int*Nteta,int*Nz,float Ea,float Eb,float Eb2)
{      
  // FILE *F;
  int i,j,g,k,save;	                	//counters
  float decal_teta;                            // gap along teta.
  float z=-L/2;
  float ez=Ea,iz=Eb ,iz2=Eb2,eo=1.85/9*2*Pi,io=Pi/2-eo;	// size of electrode contacts
  int nteta,nl,c,nez,niz,niz2,neo,nio;		//number of points
  float dteta,dl,alpha,o,dez,diz,diz2,deo,dio;
  int np=0, nt=0;				//number of points and triangles to create

  decoupe(dt,ez,&nez,&dez);
  decoupe(dt,iz,&niz,&diz);
decoupe(dt,iz2,&niz2,&diz2);
  decoupe(dt/R,eo,&neo,&deo);
  decoupe(dt/R,io,&nio,&dio);

  nteta=2*(neo+nio);				//number of points on an outer circle
  dteta=2*Pi/nteta;
  dt = dteta*R;
  decoupe(dt,(L-3*ez-3*iz - iz2)/2,&nl,&dl);
  decal_teta= -deo*neo/2;
  *Nteta=4*(neo+nio);
  *Nz=3*nez +3*niz +niz2+1;
 
  c=(int)(Pi/2*R/dt+0.5)+1;             	 //number of inner circles
  c=(c<2)?2:c;
  alpha=Pi/2/(c-1);		                //distance between inner circles
  int N[c];		                        //number of points on inner circles
  int num[c];            		       //index of first point on inner circle

  float ang[c];

  int max = (int) floor(nteta*(L/dt+2*c)*4);
    Vect3* P = new Vect3[max];
    Triangle* T = new Triangle[2*max];

  // indices and coordinates of points before electrode
  for (j=0;j<nl;j++)
  { for (i=0;i<nteta;i++)
    {
      o=i*dteta + (j-nl+1)*dteta/2 + decal_teta;
      P[np].x()=R*cos(o);
      P[np].y()=R*sin(o);
      P[np].z()=j*dl +z;
      np ++;
    }
  }
 
  // indices and coordinates of electrode points
  float zz=nl*dl +z -diz;
  for (j=0;j<=6;j++)
  {
    int gmax=(j%2==0)?((j == 4)?niz2:niz):(nez);
    if (j==0) gmax++;
    for (g=0;g<gmax;g++)
    {
      zz += (j%2==0)?((j == 4)?diz2:diz):(dez);
      for (k=0;k<4;k++)
      { for (i=0;i<neo;i++)
	{
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
  for (j=nl;j<(2*nl);j++)
  { for (i=0;i<nteta;i++)
    {
      o=i*dteta + (j-nl)*dteta/2 + decal_teta;
      P[np].x()=R*cos(o);
      P[np].y()=R*sin(o);
      P[np].z()=(j+1)*dl +z +3*ez+3*iz+iz2;
      np ++;
    }
  }


  // definition of triangles before electrode
  for (j=0;j<nl-1;j++)
  { for (i=0;i<nteta;i++)
    {
      T[nt][0]=i+nteta*j;
      T[nt][1]=(i+1)%nteta +nteta*j;
      T[nt][2]= i+nteta*(j+1);

      T[nt+1][0]=i+nteta*j;
      T[nt+1][1]=i+nteta*(j+1);
      T[nt+1][2]=(nteta+i-1)%nteta +nteta*(j+1);

      nt += 2;
    }
  }

  // definition of transition triangles
    for (i=0;i<nteta;i++)
    {
      T[nt][0]=i+nteta*(nl-1);
      T[nt][1]=(2*i+1) + nteta*nl;
      T[nt][2]=2*i + nteta*nl;
      nt++;
      T[nt][0]=i+nteta*(nl-1);
      T[nt][1]=(i+1)%nteta+nteta*(nl-1);
      T[nt][2]=(2*i+1) + nteta*nl;
      nt++;
      T[nt][0]=(i+1)%nteta+nteta*(nl-1);
      T[nt][1]=(2*i+2)%(2*nteta) + nteta*nl;
      T[nt][2]=(2*i+1) + nteta*nl;
      nt++;
    }

  // definition of electrode triangles
  for (j=0;j<(3*nez+3*niz+niz2);j++)
  {
    for (i=0;i<2*nteta;i++)
    {
      T[nt][0]=i+2*nteta*j +nteta*nl;
      T[nt][1]=i+2*nteta*(j+1) +nteta*nl;
      T[nt][2]=(2*nteta+i-1)%(2*nteta) +2*nteta*(j+1) +nteta*nl;

      T[nt+1][0]=i+2*nteta*j +nteta*nl;
      T[nt+1][1]=(i+1)%(2*nteta) +2*nteta*j +nteta*nl;
      T[nt+1][2]= i+2*nteta*(j+1) +nteta*nl;

      nt += 2;
    }
  }

  // definition of transition triangles
    for (i=0;i<nteta;i++)
    {
      T[nt][0]=i +save;
      T[nt][1]=2*i  +save-2*nteta;
      T[nt][2]=(2*i+1) +save-2*nteta;
      nt++;
      T[nt][0]=i +save;
      T[nt][1]=(2*i+1) +save-2*nteta;
      T[nt][2]=(i+1)%nteta +save;
      nt++;
      T[nt][0]=(i+1)%nteta +save;
      T[nt][1]=(2*i+1) +save-2*nteta;
      T[nt][2]=(2*i+2)%(2*nteta) +save-2*nteta;
      nt++;
    }

  // definition of triangles after electrode
  for (j=0;j<nl-1;j++)
  { for (i=0;i<nteta;i++)
    {
      T[nt][0]=i+nteta*j +save;
      T[nt][1]=(i+1)%nteta +nteta*j +save;
      T[nt][2]= i+nteta*(j+1) +save;

      T[nt+1][0]=i+nteta*j +save;
      T[nt+1][1]=i+nteta*(j+1) +save;
      T[nt+1][2]=(nteta+i-1)%nteta +nteta*(j+1) +save;

      nt += 2;
    }
  }

  save = np-nteta;
  // the lids:
  for (g=0;g<=1;g++) {
    num[c-1]=((g==0)?(0):(save));
  for (j=0;j<c;j++)
  {
    if (j==0)
	N[0]=1;
    else
	decoupe(dt/(R*sin(alpha*j)),2*Pi,&N[j],&ang[j]);
    if (j==c-1) break;
    num[j]=np;
    for (i=0;i<N[j];i++)
    {
      o=i*ang[j] + ((g==0)?((-nl+1)*dteta/2):((nl-1)*dteta/2)) + decal_teta;
      P[np].x()=R*sin(alpha*j)*cos(o);
      P[np].y()=R*sin(alpha*j)*sin(o);
      P[np].z()=g*L +z + R*cos(alpha*j) *(g*2-1) *0.3;
      np ++;
    }
  }
  for (j=2;j<c;j++)
  {
    k=0;
      T[nt][0]=num[j];
      T[nt][1]=num[j]+1;
      T[nt][2]=num[j-1];
      if (g==0) permute(&T[nt]);
      nt ++;
    for (i=1;i<N[j];i++)
    {
      if (ang[j-1]*(k+0.5) < ang[j]*(i+0.5))
      {
	T[nt][0]=num[j]+i;
	T[nt][2]=num[j-1]+k;
	T[nt][1]=num[j-1]+(k+1)%N[j-1];
	if (g==0) permute(&T[nt]);
	k=(k+1)%N[j-1];
	nt ++;
      }
	T[nt][0]=num[j]+i%N[j];
	T[nt][1]=num[j]+(i+1)%N[j];
	T[nt][2]=num[j-1]+k;
	if (g==0) permute(&T[nt]);
	nt ++;
    }
  }
  for (i=0;i<N[1];i++)
  {
    T[nt][0]=num[1]+i;
    T[nt][1]=num[1]+(i+1)%N[1];
    T[nt][2]=num[0];
    if (g==0) permute(&T[nt]);
    nt ++;
  }
  }
  Mesh surf(np,nt);
  for  (i=0;i<np;i++){
  surf.setPt(i,P[i]);
     }
  for (i=0;i<nt;i++){
    surf.setTrg(i,T[i]);
      }
  surf.save(namesurf);

  //the electrode
  if (E[0]>=0){
    // float aire= 1; //deo*R *dez /2;
    int trianglecounter=0;
    sparse_matrice stimelec(nt,12);
    for (i=0;i<(nl-1)*2*nteta+3*nteta+1;i++){
      trianglecounter++;
    }
    // int eleccounter = 0;
    for (j=0;j<=6;j++){
  	for (g=0;g<((j%2==0)?((j == 4)?niz2:niz):(nez));g++){
	 for (k=0;k<4;k++){
	  for (i=0;i<neo;i++){
	    // find for each electrode to which triangles it corresponds
	    if(j%2==1){
	      stimelec(trianglecounter,4*(j/2)+k) = 1.;
	      trianglecounter++    ;
	      stimelec(trianglecounter,4*(j/2)+k) = 1.;
	      trianglecounter++    ;
	    }
	    else{
	      trianglecounter++    ;
	      trianglecounter++    ;
	    }
	  }	      
	  for (i=0;i<nio;i++)
	    {
	      trianglecounter++;
	      trianglecounter++;
	    }
	}
	 //	if(j%2==1) eleccounter++;
	}
    }
    for (i= (nl-1)*2*nteta +3*nteta + 2*2*nteta*(3*nez+3*niz+niz2)+1;i<nt;i++){
      trianglecounter++;
    }
    stimelec.saveBin(namestim);
  }
  return 0;
}

int main(int argc, char** argv)
{
  command_usage("Make nerve geometry from existing parameters or make nerve geometry and parameter file from commandline user interface.");
    if ((!strcmp(argv[1],"-h")) | (!strcmp(argv[1],"--help"))) getHelp(argv);
    disp_argv(argc,argv);
    FILE *F,*Fgeom,*Fcond;
    // char nom[80],s[80];
    int i,Nc=2,Elec[4],E[1]={-1};
    float Ea,Eb,Eb2;
    int Nteta[Nc],Nz[Nc];
    float L[Nc],R[Nc],dt[Nc],sig[Nc+1];
    string name;
    sig[Nc]=0;
    
    if (!strcmp(argv[1], "-makeparameters")){
      // fprintf(stdout,"number of cylinders (for the moment fixed to 2) :");fflush(stdout);
      // fscanf(stdin,"%d",&Nc);
      fprintf(stdout,"size of electrode contact along z :");fflush(stdout);
      fscanf(stdin,"%f",&Ea);
      fprintf(stdout,"distance between first anode and cathode:");fflush(stdout);
      fscanf(stdin,"%f",&Eb);
      fprintf(stdout,"distance between second anode and cathode :");fflush(stdout);
      fscanf(stdin,"%f",&Eb2);

      for (i=0;i<Nc;i++)
	{  
	  if (i==0)
	    {printf(" nerve geometry (inner cylinder)\n");}
	  if (i==1)
	    {printf(" nerve geometry (outer cylinder containing electrode)\n");}
	  fprintf(stdout,"length :");fflush(stdout);
	  fscanf(stdin,"%f",&L[i]);
	  fprintf(stdout,"radius :");fflush(stdout);
	  fscanf(stdin,"%f",&R[i]);
	  fprintf(stdout,"discretization step :");fflush(stdout);
	  fscanf(stdin,"%f",&dt[i]);
	  fprintf(stdout,"conductivity :");fflush(stdout);
	  fscanf(stdin,"%f",&sig[i]);
	}
      F = fopen (argv[2],"w");
      fprintf(F,"%f\n",Ea);
      fprintf(F,"%f\n",Eb);
      fprintf(F,"%f\n",Eb2);
      for (i=0;i<Nc;i++)
	{  
	  fprintf(F,"%f\n",L[i]);
	  fprintf(F,"%f\n",R[i]);
	  fprintf(F,"%f\n",dt[i]);
	  fprintf(F,"%f\n",sig[i]);	
	}  
      fclose(F);
    }
    else if (!strcmp(argv[1],"-useparameters")){
      // Read parameter file :
      F = fopen (argv[2],"r");
      if (F==NULL) {perror ("erreur fopen: \n"); return -1;}
      fscanf(F,"%f",&Ea);
      fscanf(F,"%f",&Eb);;
      fscanf(F,"%f",&Eb2);
      for (i=0;i<Nc;i++)
	{  
	  fscanf(F,"%f",&L[i]);
	  fscanf(F,"%f",&R[i]);
	  fscanf(F,"%f",&dt[i]);
	  fscanf(F,"%f",&sig[i]);
	  /*	  fclose(F);*/
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
    Fgeom = fopen (argv[3],"w");
    Fcond = fopen (argv[4],"w");
    if (Fgeom==NULL) {perror ("erreur fopen: .geom\n"); return -1;}
    if (Fcond==NULL) {perror ("erreur fopen: .cond\n"); return -1;}
    fprintf(Fgeom,"# Domain Description 1.0\n\n");
    fprintf(Fgeom,"Interfaces %d Mesh\n\n",Nc);
    fprintf(Fcond,"# Properties Description 1.0 (Conductivities)\n\n");
    for (i=0;i<Nc;i++){ 
      cylindre(argv[5+i],argv[7],L[i],R[i],dt[i],(i==Nc-1)?(Elec):(E),&Nteta[i],&Nz[i],Ea,Eb,Eb2);
      name = argv[5+i];
      name = name.substr(name.rfind("/")+1); // only keep the file name without the directories (after the last /)
      fprintf(Fgeom,"%s\n",name.c_str());    // c_str to convert string back to char
      }
    fprintf(Fcond,"Air        0.0\n");
    fprintf(Fcond,"CSF        %f\n",sig[1]);
    fprintf(Fcond,"Nerve      %f\n",sig[0]);
    fprintf(Fgeom,"\nDomains %d\n\n",Nc+1);
    fprintf(Fgeom,"Domain CSF 1 -2\n");
    fprintf(Fgeom,"Domain Nerve -1 shared\n");
    fprintf(Fgeom,"Domain Air 2\n");
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
    cout << "               Output stimelec sparse matrix" << endl << endl;
    cout << "         -makeparameters       " <<endl;    
    cout << "               Output parameter file" << endl;
    cout << "               Output geom file " << endl;
    cout << "               Output cond file " << endl;
    cout << "               Output surf1 tri mesh " << endl;
    cout << "               Output surf2 tri mesh " << endl;
    cout << "               Output stimelec sparse matrix" << endl << endl;

    exit(0);
}
