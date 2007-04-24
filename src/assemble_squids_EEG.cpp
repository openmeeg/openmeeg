#include "mainHeader.h"
#include "danielsson.h"
#include "operateurs.h"

int* computeVindexes(geometry &geo,int* n_indexes=0)
{
    int count=0;
    for(int i=0;i<geo.nb();i++) count+=geo.getM(i).nbr_pts();

    int *ret=new int[count];
    if(n_indexes!=0) *n_indexes=count;

    count=0;
    int offset=0;
    for(int i=0;i<geo.nb();i++)
    {
        for(int j=0;j<geo.getM(i).nbr_pts();j++) {ret[count]=count+offset; count++;}
        offset+=geo.getM(i).nbr_trg();
    }

    return ret;
}

void assemble_xToEEGresponse( geometry &geo, matrice &mat, const matrice &positions )
//EEG patches positions are reported line by line in the positions matrix
//mat is supposed to be filled with zeros
//mat is the linear application which maps x (the unknown vector in symmetric system) -> v (potential at the electrodes)
{
    mesh &extLayer=geo.getM(geo.nb()-1);

    //first we calculate the offset of the potential on the external layer in the unknown vector x
    int offset=0;
    for(int l=0;l<geo.nb()-1;l++)
    {
        offset+=geo.getM(l).nbr_pts(); offset+=geo.getM(l).nbr_trg();
    }

    vect3 current_position;
    vect3 current_alphas;
    int current_nearestNumber;
    for(size_t i=0;i<positions.nlin();i++)
    {
        for(int k=0;k<3;k++) current_position[k]=positions(i,k);
        dist_point_mesh(current_position ,extLayer,current_alphas,current_nearestNumber);
        mat(i,extLayer.trngl(current_nearestNumber).id1()+offset)=current_alphas[0];
        mat(i,extLayer.trngl(current_nearestNumber).id2()+offset)=current_alphas[1];
        mat(i,extLayer.trngl(current_nearestNumber).id3()+offset)=current_alphas[2];
    }

}

void assemble_xToMEGresponseContrib( geometry &geo, matrice &mat, const matrice &positions, const matrice &orientations )
//MEG patches positions are reported line by line in the positions matrix (same for positions)
//mat is supposed to be filled with zeros
//mat is the linear application which maps x (the unknown vector in symmetric system) -> bFerguson (contrib to MEG response)
{
    matrice myFergusonMatrix(3*mat.nlin(),mat.ncol());
    myFergusonMatrix.set(0.0);
    const int nsquids=(int)positions.nlin();
    vect3 *positionsVectArray=new vect3[nsquids];

    int n_indexes;
    int* vIndexes=computeVindexes(geo,&n_indexes);

    for(int i=0;i<nsquids;i++)
    {
        positionsVectArray[i][0]=positions(i,0);
        positionsVectArray[i][1]=positions(i,1);
        positionsVectArray[i][2]=positions(i,2);
    }

    assemble_ferguson(geo,myFergusonMatrix,positionsVectArray,nsquids);

    for(size_t i=0;i<mat.nlin();i++)
    {
        for(int j=0;j<n_indexes;j++)
        {
            vect3 fergusonField(myFergusonMatrix(3*i,vIndexes[j]),myFergusonMatrix(3*i+1,vIndexes[j]),myFergusonMatrix(3*i+2,vIndexes[j]));
            vect3 normalizedDirection(orientations(i,0),orientations(i,1),orientations(i,2));
            normalizedDirection.normalize();
            mat(i,vIndexes[j])=fergusonField*normalizedDirection;
        }
    }

    delete[] positionsVectArray;
    delete[] vIndexes;
}


void assemble_sToMEGresponseContrib( mesh &sources_mesh, matrice &mat, const matrice &positions, const matrice &orientations )
//MEG patches positions are reported line by line in the positions matrix (same for positions)
//mat is supposed to be filled with zeros
//mat is the linear application which maps x (the unknown vector in symmetric system) -> binf (contrib to MEG response)
{
    matrice myFergusonMatrix(3*mat.nlin(),mat.ncol());
    myFergusonMatrix.set(0.0);
    const int nsquids=(int)positions.nlin();
    vect3 *positionsVectArray=new vect3[nsquids];

    for(int i=0;i<nsquids;i++)
    {
        positionsVectArray[i][0]=positions(i,0);
        positionsVectArray[i][1]=positions(i,1);
        positionsVectArray[i][2]=positions(i,2);
    }

    for(size_t i=0;i<mat.nlin();i++) operateurFerguson(positionsVectArray[i],sources_mesh,myFergusonMatrix,3*(int)i,0);

    for(size_t i=0;i<mat.nlin();i++)
    {
        for(size_t j=0;j<mat.ncol();j++)
        {
            vect3 fergusonField(myFergusonMatrix(3*i,j),myFergusonMatrix(3*i+1,j),myFergusonMatrix(3*i+2,j));
            vect3 normalizedDirection(orientations(i,0),orientations(i,1),orientations(i,2));
            normalizedDirection.normalize();
            mat(i,j)=fergusonField*normalizedDirection;
        }
    }

    delete[] positionsVectArray;
}


// creates the S2MEG matrix with unconstrained orientations for the sources.
void assemble_sToMEGresponseContrib_point( matrice&dipoles, matrice &mat, const matrice &positions, const matrice &orientations )
//MEG patches positions are reported line by line in the positions matrix (same for positions)
//mat is supposed to be filled with zeros
//mat is the linear application which maps x (the unknown vector in symmetric system) -> binf (contrib to MEG response)
//sources is the name of a file containing the description of the sources - one dipole per line: x1 x2 x3 n1 n2 n3, x being the position and n the orientation.
{
    if(dipoles.ncol()!=6) {std::cerr<<"Dipoles File Format Error"<<std::endl; exit(1);}
    int nd=(int)dipoles.nlin();
    std::vector<vect3> Rs,Qs;
    for(int i=0;i<nd;i++)
    {
        vect3 r(3),q(3);
        for(int j=0;j<3;j++) r[j]=dipoles(i,j);
        for(int j=3;j<6;j++) q[j-3]=dipoles(i,j);
        Rs.push_back(r); Qs.push_back(q);
    }
  //Rs and Qs respectively contains positions and orientations of the dipoles.

  //this matrix will contain the field generated at the location of the i-th squid by the j-th source
    matrice SignalMatrix(3*mat.nlin(),mat.ncol());
    SignalMatrix.set(0.0);
    const int nsquids=(int)positions.nlin();
    vect3 *positionsVectArray=new vect3[nsquids];

    for(int i=0;i<nsquids;i++)
    {
        positionsVectArray[i][0]=positions(i,0);
        positionsVectArray[i][1]=positions(i,1);
        positionsVectArray[i][2]=positions(i,2);
    }

  // the following routine is the equivalent of operateurFerguson for pointlike dipoles.
    for(size_t i=0;i<mat.nlin();i++)
    {
        for(unsigned int j=0;j<0+mat.ncol();j++)
        {
            vect3 diff=positionsVectArray[i]-Rs[j];
            double norm_diff=diff.norme();
            vect3 v=diff/(norm_diff*norm_diff*norm_diff)^Qs[j];

            SignalMatrix(3*i+0,j)=v.X();
            SignalMatrix(3*i+1,j)=v.Y();
            SignalMatrix(3*i+2,j)=v.Z();
        }
    }


    for(size_t i=0;i<mat.nlin();i++)
    {
        for(size_t j=0;j<mat.ncol();j++)
        {
            vect3 fergusonField(SignalMatrix(3*i,j),SignalMatrix(3*i+1,j),SignalMatrix(3*i+2,j));
            vect3 normalizedDirection(orientations(i,0),orientations(i,1),orientations(i,2));
            normalizedDirection.normalize();
            mat(i,j)=fergusonField*normalizedDirection/(4*M_PI);
        }
    }


    delete[] positionsVectArray;
}

