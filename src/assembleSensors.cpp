/* FILE: $Id$ */

/*
Project Name : OpenMEEG

author            : $Author$
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

#include "assemble.h"
#include "danielsson.h"
#include "operators.h"
#include "sensors.h"

void assemble_ferguson(const Geometry &geo, matrice &mat, const Vect3* pts,const int n);

int* computeVindexes(const Geometry &geo,int* n_indexes=0)
{
    int count=0;
    for(int i=0;i<geo.nb();i++) count+=geo.getM(i).nbPts();

    int *ret=new int[count];
    if(n_indexes!=0) *n_indexes=count;

    count=0;
    int offset=0;
    for(int i=0;i<geo.nb();i++)
    {
        for(int j=0;j<geo.getM(i).nbPts();j++) {ret[count]=count+offset; count++;}
        offset+=geo.getM(i).nbTrgs();
    }

    return ret;
}

void assemble_vToEEG(matrice &mat, const Geometry &geo, const matrice &positions )
//EEG patches positions are reported line by line in the positions matrix
//mat is supposed to be filled with zeros
//mat is the linear application which maps x (the unknown vector in symmetric system) -> v (potential at the electrodes)
{
    int newsize = geo.size()-(geo.getM(geo.nb()-1)).nbTrgs();
    mat = matrice(positions.nlin(),newsize);
    mat.set(0.0);

    const Mesh& extLayer=geo.getM(geo.nb()-1);

    //first we calculate the offset of the potential on the external layer in the unknown vector x
    int offset=0;
    for(int l=0; l < geo.nb()-1; l++)
    {
        offset += geo.getM(l).nbPts(); offset+=geo.getM(l).nbTrgs();
    }

    Vect3 current_position;
    Vect3 current_alphas;
    int current_nearestNumber;
    for(size_t i=0; i < positions.nlin(); i++)
    {
        for(int k=0;k<3;k++) current_position(k)=positions(i,k);
        dist_point_mesh(current_position ,extLayer,current_alphas,current_nearestNumber);
        mat(i,extLayer.getTrg(current_nearestNumber).s1()+offset)=current_alphas(0);
        mat(i,extLayer.getTrg(current_nearestNumber).s2()+offset)=current_alphas(1);
        mat(i,extLayer.getTrg(current_nearestNumber).s3()+offset)=current_alphas(2);
    }

}

vToEEG_matrice::vToEEG_matrice(const Geometry &geo, const matrice &positions) {
    assemble_vToEEG(*this, geo, positions);
}

void assemble_vToMEG(matrice &mat, const Geometry &geo, const Sensors &sensors)
//MEG patches positions are reported line by line in the positions matrix (same for positions)
//mat is supposed to be filled with zeros
//mat is the linear application which maps x (the unknown vector in symmetric system) -> bFerguson (contrib to MEG response)
{
    matrice positions = sensors.getPositions();
    matrice orientations = sensors.getOrientations();

    int newsize = geo.size()-(geo.getM(geo.nb()-1)).nbTrgs();
    mat = matrice(positions.nlin(),newsize);
    mat.set(0.0);

    matrice myFergusonMatrix(3*mat.nlin(),mat.ncol());
    myFergusonMatrix.set(0.0);
    const int nsquids=(int)positions.nlin();
    Vect3 *positionsVectArray=new Vect3[nsquids];

    int n_indexes;
    int* vIndexes=computeVindexes(geo,&n_indexes);

    for(int i=0;i<nsquids;i++)
    {
        positionsVectArray[i](0) = positions(i,0);
        positionsVectArray[i](1) = positions(i,1);
        positionsVectArray[i](2) = positions(i,2);
    }

    assemble_ferguson(geo,myFergusonMatrix,positionsVectArray,nsquids);

    for(size_t i=0;i<mat.nlin();i++)
    {
        progressbar(i,mat.nlin());
        #ifdef USE_OMP
        #pragma omp parallel for
        #endif
        for(int j=0;j<n_indexes;j++)
        {
            Vect3 fergusonField(myFergusonMatrix(3*i,vIndexes[j]),myFergusonMatrix(3*i+1,vIndexes[j]),myFergusonMatrix(3*i+2,vIndexes[j]));
            Vect3 normalizedDirection(orientations(i,0),orientations(i,1),orientations(i,2));
            normalizedDirection.normalize();
            mat(i,vIndexes[j]) = fergusonField*normalizedDirection;
        }
    }

    delete[] positionsVectArray;
    delete[] vIndexes;
}

vToMEG_matrice::vToMEG_matrice(const Geometry &geo, const Sensors &sensors) {
    assemble_vToMEG(*this, geo, sensors);
}

void assemble_sToMEG(matrice &mat, const Mesh &sources_mesh, const Sensors &sensors)
//MEG patches positions are reported line by line in the positions matrix (same for positions)
//mat is supposed to be filled with zeros
//mat is the linear application which maps x (the unknown vector in symmetric system) -> binf (contrib to MEG response)
{
    matrice positions = sensors.getPositions();
    matrice orientations = sensors.getPositions();

    mat = matrice(positions.nlin(),sources_mesh.nbPts());
    mat.set(0.0);

    matrice myFergusonMatrix(3*mat.nlin(),mat.ncol());
    myFergusonMatrix.set(0.0);
    const int nsquids=(int)positions.nlin();
    Vect3 *positionsVectArray=new Vect3[nsquids];

    for(int i=0;i<nsquids;i++)
    {
        positionsVectArray[i](0)=positions(i,0);
        positionsVectArray[i](1)=positions(i,1);
        positionsVectArray[i](2)=positions(i,2);
    }

    for(size_t i=0;i<mat.nlin();i++) {
        progressbar(i,mat.nlin());
        operatorFerguson(positionsVectArray[i],sources_mesh,myFergusonMatrix,3*(int)i,0);
    }
    
    for(size_t i=0;i<mat.nlin();i++)
    {
        for(size_t j=0;j<mat.ncol();j++)
        {
            Vect3 fergusonField(myFergusonMatrix(3*i,j),myFergusonMatrix(3*i+1,j),myFergusonMatrix(3*i+2,j));
            Vect3 normalizedDirection(orientations(i,0),orientations(i,1),orientations(i,2));
            normalizedDirection.normalize();
            mat(i,j)=fergusonField*normalizedDirection;
        }
    }

    delete[] positionsVectArray;
}

sToMEG_matrice::sToMEG_matrice(const Mesh &sources_mesh, const Sensors &sensors) {
    assemble_sToMEG(*this, sources_mesh, sensors);
}

// creates the S2MEG matrix with unconstrained orientations for the sources.
void assemble_sToMEGdip( matrice &mat, const matrice& dipoles, const Sensors &sensors)
//MEG patches positions are reported line by line in the positions matrix (same for positions)
//mat is supposed to be filled with zeros
//mat is the linear application which maps x (the unknown vector in symmetric system) -> binf (contrib to MEG response)
//sources is the name of a file containing the description of the sources - one dipole per line: x1 x2 x3 n1 n2 n3, x being the position and n the orientation.
{
    matrice positions = sensors.getPositions();
    matrice orientations = sensors.getOrientations();

    mat = matrice(positions.nlin(),dipoles.nlin());

    if(dipoles.ncol()!=6) {std::cerr<<"Dipoles File Format Error"<<std::endl; exit(1);}
    int nd=(int)dipoles.nlin();
    std::vector<Vect3> Rs,Qs;
    for(int i=0;i<nd;i++)
    {
        Vect3 r(3),q(3);
        for(int j=0;j<3;j++) r(j)=dipoles(i,j);
        for(int j=3;j<6;j++) q(j-3)=dipoles(i,j);
        Rs.push_back(r); Qs.push_back(q);
    }
    //Rs and Qs respectively contains positions and orientations of the dipoles.

    //this matrix will contain the field generated at the location of the i-th squid by the j-th source
    matrice SignalMatrix(3*mat.nlin(),mat.ncol());
    SignalMatrix.set(0.0);
    const int nsquids=(int)positions.nlin();
    Vect3 *positionsVectArray=new Vect3[nsquids];

    for(int i=0;i<nsquids;i++)
    {
        positionsVectArray[i](0)=positions(i,0);
        positionsVectArray[i](1)=positions(i,1);
        positionsVectArray[i](2)=positions(i,2);
    }

    // the following routine is the equivalent of operatorFerguson for pointlike dipoles.
    for(size_t i=0;i<mat.nlin();i++)
    {
        for(unsigned int j=0;j<0+mat.ncol();j++)
        {
            Vect3 diff=positionsVectArray[i]-Rs[j];
            double norm_diff=diff.norme();
            Vect3 v = Qs[j] ^ diff/(norm_diff*norm_diff*norm_diff);

            SignalMatrix(3*i+0,j) = v.x();
            SignalMatrix(3*i+1,j) = v.y();
            SignalMatrix(3*i+2,j) = v.z();
        }
    }

    for(size_t i=0;i<mat.nlin();i++)
    {
        for(size_t j=0;j<mat.ncol();j++)
        {
            Vect3 fergusonField(SignalMatrix(3*i,j),SignalMatrix(3*i+1,j),SignalMatrix(3*i+2,j));
            Vect3 normalizedDirection(orientations(i,0),orientations(i,1),orientations(i,2));
            normalizedDirection.normalize();
            mat(i,j)=fergusonField*normalizedDirection/(4*M_PI);
        }
    }

    delete[] positionsVectArray;
}

sToMEGdip_matrice::sToMEGdip_matrice(const matrice &dipoles, const Sensors &sensors) {
    assemble_sToMEGdip(*this, dipoles, sensors);
}


