/* FILE: $Id: assembleSensors.cpp 222 2008-04-08 06:14:41Z gramfort $ */

/*
Project Name : OpenMEEG

version           : $Revision: 222 $
last revision     : $Date: 2008-04-08 08:14:41 +0200 (Tue, 08 Apr 2008) $
modified by       : $LastChangedBy: gramfort $
last modified     : $LastChangedDate: 2008-04-08 08:14:41 +0200 (Tue, 08 Apr 2008) $

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

#include "assemble.h"
#include "danielsson.h"
#include "operators.h"
#include "sensors.h"

void assemble_ferguson(const Geometry &geo, Matrix &mat, const Vect3* pts,const int n);

// EEG patches positions are reported line by line in the positions Matrix
// mat is supposed to be filled with zeros
// mat is the linear application which maps x (the unknown vector in symmetric system) -> v (potential at the electrodes)
void assemble_Head2EEG(SparseMatrix &mat, const Geometry &geo, const Matrix &positions )
{
    int newsize = geo.size()-(geo.getM(geo.nb()-1)).nbTrgs();
    mat = SparseMatrix(positions.nlin(),newsize);

    const Mesh& extLayer = geo.getM(geo.nb()-1);

    //first we calculate the offset of the potential on the external layer in the unknown vector x
    int offset=0;
    for(int l=0; l < geo.nb()-1; l++)
    {
        offset += geo.getM(l).nbPts();
        offset += geo.getM(l).nbTrgs();
    }

    Vect3 current_position;
    Vect3 current_alphas;
    int current_nearestNumber;
    for(size_t i=0; i < positions.nlin(); i++)
    {
        for(int k=0;k<3;k++) current_position(k) = positions(i,k);
        dist_point_mesh(current_position, extLayer, current_alphas, current_nearestNumber);
        mat(i,extLayer.getTrg(current_nearestNumber).s1()+offset) = current_alphas(0);
        mat(i,extLayer.getTrg(current_nearestNumber).s2()+offset) = current_alphas(1);
        mat(i,extLayer.getTrg(current_nearestNumber).s3()+offset) = current_alphas(2);
    }
}

Head2EEGMat::Head2EEGMat(const Geometry &geo, const Matrix &positions) {
    assemble_Head2EEG(*this, geo, positions);
}

// MEG patches positions are reported line by line in the positions Matrix (same for positions)
// mat is supposed to be filled with zeros
// mat is the linear application which maps x (the unknown vector in symmetric system) -> bFerguson (contrib to MEG response)
void assemble_Head2MEG(Matrix &mat, const Geometry &geo, const Sensors &sensors)
{
    Matrix positions = sensors.getPositions();
    Matrix orientations = sensors.getOrientations();
    const size_t nbIntegrationPoints = sensors.getNumberOfPositions();
    int p0_p1_size = geo.size()-(geo.getM(geo.nb()-1)).nbTrgs();
    int geo_number_points = geo.getNumberOfPoints();

    mat = Matrix(nbIntegrationPoints,p0_p1_size);
    mat.set(0.0);

    Matrix myFergusonMatrix(3*nbIntegrationPoints,geo_number_points);
    myFergusonMatrix.set(0.0);

    Vect3 *positionsVectArray = new Vect3[nbIntegrationPoints];
    for(size_t i=0;i<nbIntegrationPoints;i++)
    {
        positionsVectArray[i](0) = positions(i,0);
        positionsVectArray[i](1) = positions(i,1);
        positionsVectArray[i](2) = positions(i,2);
    }

    assemble_ferguson(geo,myFergusonMatrix,positionsVectArray,nbIntegrationPoints);

    // Compute indexes of V indexes (P1 elements)
    int* vindex = new int[geo_number_points];
    int count = 0;
    int offset = 0;
    for(int i=0;i<geo.nb();i++)
    {
        for(int j=0;j<geo.getM(i).nbPts();j++) {
            vindex[count] = count+offset;
            count++;
        }
        offset += geo.getM(i).nbTrgs();
    }

    for(size_t i=0;i<nbIntegrationPoints;i++)
    {
        progressbar(i,nbIntegrationPoints);
        #ifdef USE_OMP
        #pragma omp parallel for
        #endif
        for(int j=0;j<geo_number_points;j++)
        {
            Vect3 fergusonField(myFergusonMatrix(3*i,j),myFergusonMatrix(3*i+1,j),myFergusonMatrix(3*i+2,j));
            Vect3 normalizedDirection(orientations(i,0),orientations(i,1),orientations(i,2));
            normalizedDirection.normalize();
            mat(i,vindex[j]) = fergusonField*normalizedDirection;
        }
    }

    mat = sensors.getWeightsMatrix() * mat; // Apply weights

    delete[] positionsVectArray;
    delete[] vindex;
}

Head2MEGMat::Head2MEGMat(const Geometry &geo, const Sensors &sensors) {
    assemble_Head2MEG(*this, geo, sensors);
}

// MEG patches positions are reported line by line in the positions Matrix (same for positions)
// mat is supposed to be filled with zeros
// mat is the linear application which maps x (the unknown vector in symmetric system) -> binf (contrib to MEG response)
void assemble_SurfSource2MEG(Matrix &mat, const Mesh &sources_mesh, const Sensors &sensors)
{
    Matrix positions = sensors.getPositions();
    Matrix orientations = sensors.getOrientations();

    mat = Matrix(positions.nlin(),sources_mesh.nbPts());
    mat.set(0.0);

    Matrix myFergusonMatrix(3*mat.nlin(),mat.ncol());
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

    mat = sensors.getWeightsMatrix() * mat; // Apply weights

    delete[] positionsVectArray;
}

SurfSource2MEGMat::SurfSource2MEGMat(const Mesh &sources_mesh, const Sensors &sensors) {
    assemble_SurfSource2MEG(*this, sources_mesh, sensors);
}

// creates the DipSource2MEG Matrix with unconstrained orientations for the sources.
//MEG patches positions are reported line by line in the positions Matrix (same for positions)
//mat is supposed to be filled with zeros
//sources is the name of a file containing the description of the sources - one dipole per line: x1 x2 x3 n1 n2 n3, x being the position and n the orientation.
void assemble_DipSource2MEG( Matrix &mat, const Matrix& dipoles, const Sensors &sensors)
{
    Matrix positions = sensors.getPositions();
    Matrix orientations = sensors.getOrientations();

    mat = Matrix(positions.nlin(),dipoles.nlin());

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

    //this Matrix will contain the field generated at the location of the i-th squid by the j-th source
    Matrix SignalMatrix(3*mat.nlin(),mat.ncol());
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
        for(unsigned int j=0;j<mat.ncol();j++)
        {
            Vect3 diff = positionsVectArray[i]-Rs[j];
            double norm_diff = diff.norm();
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

    mat = sensors.getWeightsMatrix() * mat; // Apply weights

    delete[] positionsVectArray;
}

DipSource2MEGMat::DipSource2MEGMat(const Matrix &dipoles, const Sensors &sensors) {
    assemble_DipSource2MEG(*this, dipoles, sensors);
}


